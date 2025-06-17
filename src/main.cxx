#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <deque>
#include <climits>
#include <vector>
#include <cstdint>
#include "bits.hxx"

/*
 *
 * Graph Coloring for Register Allocation
 *
 * Algorithm explained in the book "modern compiler implementation in C"
 * 	by Andrew W. Appel
 *
 * Not a 1:1 mapping because I only followed the explanations,
 * 	not the implementation and the decisions that came with it.
 * All data structures are thus not necessarily using the best known
 * 	representation, just what made the most sense to me at conception.
 *
 * This code was intended to be plugged in a compiler project of mine
 * 	however with the time that passed I feel like I would rather make a
 * 	new one from scratch and so this would also see modifications
 * 	most likely in the API to be more flexible with lowered instructions
 * 	since we want to go from a universal bytecode to an architecture-
 * 	aware bytecode, and also eventually supporting allocations for
 * 	different types (i8-i64, f32-f64, struct, simd vectors, etc)
 *
 */

using reg = std::uint32_t;

struct instr
{
	enum opcode_t : reg {
		def,   // rd <- ?
		req,   // rd -> ? // ? <- rs1
		copy,  // rd <- rs1
		load,  // rd <- [rs1]
		store, // rd -> [rs1] // TODO: [rs1] <- rs2
		load_local,  // rd <- local_#rs1 // rd <- imm1[fp]
		store_local, // rd -> local_#rs1 // TODO: imm1[fp] <- rs2
		add,   // rd <- rs1 + rs2
		imm,   // rd <- #rs1 // rd <- imm1
		umul,  // rd <- rs1 * rs2

		clobber, // rd <- undef because of rs1
	} opcode;
	reg rd ;
	reg rs1; // imm1
	reg rs2;
};

reg lsb(reg r)
{
	return r & -r;
}

reg ctz(reg r)
{
	return __builtin_ctz(r);
}

// This function is called with a mask of 'free' registers and
// 	selects any of them. This is a customization point where
// 	a different heuristic could be used
reg some_bit_index(reg r)
{
	return ctz(lsb(r));
}

reg bit(reg r)
{
	assert(r < CHAR_BIT * sizeof r);
	return reg{1} << r;
}

size_t round_up(size_t from, size_t to)
{
	return (from + to - 1) / to * to;
}

struct graph
{
	bits::managed nonmove_;
	bits::managed move_;
	bits::managed removed;
	bits::managed frozen;

	graph(reg count);

	void nonmove(reg s, reg t);
	void move(reg s, reg t);
	void remove(reg s);
	void freeze(reg s);

	bool is_free(reg s, reg t) const;
	bool is_move(reg s, reg t) const;
	bool is_nonmove(reg s, reg t) const;
	bool is_frozen(reg s) const;
	bool is_live(reg s) const;
	bool is_nonmove(reg s) const;
	size_t idx(reg s, reg t) const;
	size_t order() const;
	size_t deg(reg s) const;
	size_t move_related(reg s) const;
};

graph::graph(reg count)
	: nonmove_((count = round_up(count, bits::managed::unit_size()), count * count)),
	  move_   (count * count),
	  removed (count),
	  frozen  (count)
{
}

size_t graph::idx(reg s, reg t) const
{
	return s * order() + t;
}

size_t graph::order() const
{
	return removed.size * bits::unit_bits;
}

void graph::nonmove(reg s, reg t)
{
	nonmove_.set(idx(s, t));
	nonmove_.set(idx(t, s));
}

void graph::move(reg s, reg t)
{
	move_.set(idx(s, t));
	move_.set(idx(t, s));
}

bool graph::is_free(reg s, reg t) const
{
	return !is_move(s, t) & !is_nonmove(s, t);
}

bool graph::is_move(reg s, reg t) const
{
	return move_[idx(s, t)] & !is_nonmove(s, t);
}

bool graph::is_nonmove(reg s, reg t) const
{
	return nonmove_[idx(s, t)];
}

void graph::remove(reg s)
{
	removed.set(s);
}

void graph::freeze(reg s)
{
	frozen.set(s);
}

size_t graph::deg(reg s) const
{
	const auto offs = idx(s, 0);
	return bits::count(
		( move_.lazy(offs) | nonmove_.lazy(offs) ) & ~removed.lazy()
	).eval(order());
}

size_t graph::move_related(reg s) const
{
	const auto offs = idx(s, 0);
	return bits::count(
		move_.lazy(offs) & ~(nonmove_.lazy(offs) | removed.lazy())
	).eval(order());
}

bool graph::is_nonmove(reg s) const
{
	return is_live(s) && (is_frozen(s) || !move_related(s));
}

bool graph::is_frozen(reg s) const
{
	return frozen[s];
}

bool graph::is_live(reg s) const
{
	return !removed[s];
}

template <typename T>
std::ostream &operator<<(std::ostream &os, std::vector<T> const &v)
{
	os << "{\n";
	for (size_t i = 0; i < v.size(); ++i) {
		os << " " << i << ":" << v[i];
	}
	os << "}\n";
	return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, std::deque<T> const &q)
{
	os << "{  ";
	for (size_t i = 0; i < q.size(); ++i) {
		os << q[i] << "  ";
	}
	os << "}\n";
	return os;
}

std::ostream &operator<<(std::ostream &os, bits::managed const &b)
{
	os << "| ";
	for (size_t i = 0; i < b.size * bits::unit_bits; ++i) {
		if (b[i]) {
			os << i << ' ';
		}
	}
	os << "| ";
	return os;
}

std::ostream &operator<<(std::ostream &os, instr const &ins)
{
	switch (ins.opcode) {
	case instr::def:
		return os << "def   R." << ins.rd << '\n';
	case instr::req:
		return os << "req   R." << ins.rd << '\n';
	case instr::copy:
		return os << "copy  R." << ins.rd << ", R." << ins.rs1 << '\n';
	case instr::load:
		return os << "load  R." << ins.rd << ", [R." << ins.rs1 << "]" << '\n';
	case instr::store:
		return os << "store R." << ins.rd << ", [R." << ins.rs1 << "]" << '\n';
	case instr::add:
		return os << "add   R." << ins.rd << ", R." << ins.rs1 << ", R." << ins.rs2 << '\n';
	case instr::umul:
		return os << "umul  R." << ins.rd << ", R." << ins.rs1 << ", R." << ins.rs2 << '\n';
	case instr::imm:
		return os << "imm   R." << ins.rd << ", #" << ins.rs1 << '\n';
	case instr::load_local:
		return os << "load  R." << ins.rd << ", L." << ins.rs1 << "[fp]" << '\n';
	case instr::store_local:
		return os << "store R." << ins.rd << ", L." << ins.rs1 << "[fp]" << '\n';
	case instr::clobber:
		return os << "// R." << ins.rd << " clobbered by R." << ins.rs1 << '\n';
	default:
		assert(false);
	}
}

std::ostream &operator<<(std::ostream &os, graph const &g)
{
	for (reg r = 0; r < g.order(); ++r) {
		if (!g.is_live(r))
			continue;
		os << (g.is_frozen(r) ? 'F': ' ');
		os << r << "(d=" << g.deg(r) << ",m=" << g.move_related(r) << "): ";
		bool first = true;
		for (reg t = 0; t < g.order(); ++t) {
			if (!g.is_live(t) || g.is_free(r, t))
				continue;
			if (first) {
				first = false;
			} else {
				os << " ";
			}
			os << (g.is_move(r, t) ? 'M': ' ') << t;
		}
		os << '\n';
	}
	os << '\n';
	return os;
}

template <typename T>
class stack : public std::deque<T>
{
public:
	void push(T val);
	T pop();
};

template <typename T>
void stack<T>::push(T val)
{
	this->push_back(val);
}

template <typename T>
T stack<T>::pop()
{
	T val = this->back();
	this->pop_back();
	return val;
}

struct code_t : public std::vector<instr>
{
	reg phys_regs;
	reg virt_regs;

	reg regs() const { return phys_regs + virt_regs; }
};

code_t canonical(code_t const &code)
{
	code_t canon{ {}, code.phys_regs, code.virt_regs };
	for (auto ins : code) {
		switch (ins.opcode) {
		case instr::def:
		case instr::req:
		case instr::copy:
		case instr::load:
		case instr::store:
		case instr::load_local:
		case instr::store_local:
		case instr::imm:
		case instr::clobber:
			canon.push_back(ins);
			break;
		case instr::add:
			canon.emplace_back(instr{ instr::copy, ins.rd, ins.rs1 });
			canon.emplace_back(instr{ instr::add , ins.rd, ins.rd , ins.rs2 });
			break;
		case instr::umul:
			canon.emplace_back(instr{ instr::copy, 2, ins.rs1 });
			canon.emplace_back(instr{ instr::umul, 2, 2, ins.rs2 });
			canon.emplace_back(instr{ instr::clobber, 0, ins.rd });
			canon.emplace_back(instr{ instr::copy, ins.rd, 2 });
			break;
		}
	}
	return canon;
}

graph gen_graph(code_t const &code)
{
	// first registers correspond to the physical registers and all interfere with each other
	// when strictly necessary, instructions can address physical registers
	// this lets you allocate once you know the target platform
	// for instance on x86-64 unsigned multiply is MUL %eax, %eax, r/m32 [clobbers %edx]
	// this also useful for calling conventions,
	// if parameters are r0, r1, r2, and return value in r3, r4,
	// you can allocate them seamlessly with
	// copy r0, %edi; copy r1, %esi; copy r2, %edx;
	// copy %eax, r3; copy %edx, r4;
	const auto regs = code.regs();
	const auto max_index = round_up(1+code.size(), bits::unit_bits);
	graph g(regs);
	// [definition, last use)
	bits::managed live(regs * max_index);
	const auto idx = [&] (reg r, size_t index) { return r * max_index + index; };
	const auto use = [&] (reg r, size_t index) { live.set  (idx(r, index)); };
	const auto def = [&] (reg r, size_t index) { live.clear(idx(r, index)); };
	const auto clobber = [&] (reg virt, reg phys) { g.nonmove(virt, phys); };

	static int _iter = 0;
	std::cout << "iter #" << _iter++ << ":\n";

	size_t i = code.size() - 1;
	do {
		for (reg r = 0; r < regs; ++r) {
			live.assign(idx(r, i), live[ idx(r, i+1) ]);
		}

		const auto ins = code[i];
		switch (ins.opcode) {
		case instr::copy:
			def(ins.rd , i);
			use(ins.rs1, i);
			g.move(ins.rd, ins.rs1);
			break;
		case instr::imm:
		case instr::load_local:
			def(ins.rd, i);
			break;
		case instr::load:
			def(ins.rd, i);
			use(ins.rs1, i);
			break;
		case instr::add:
		case instr::umul:
			def(ins.rd, i);
			use(ins.rs1, i);
			use(ins.rs2, i);
			break;
		case instr::store:
			use(ins.rd, i);
			use(ins.rs1, i);
			break;
		case instr::store_local:
			use(ins.rd, i);
			break;
		case instr::def:
			def(ins.rd, i);
			break;
		case instr::req:
			use(ins.rd, i);
			break;
		case instr::clobber:
			clobber(ins.rd, ins.rs1);
			break;
		}
	} while (i--);

	for (reg i = 0; i < code.size(); ++i) {
		bits::managed _live(regs);
		for (reg r = 0; r < regs; ++r) {
			if (live[idx(r, i)])
				_live.set(r);
		}
		std::cout << i << ": " << _live << code[i];
	}

	// render physical registers
	for (reg t = 1; t < code.phys_regs; ++t) {
		for (reg s = 0; s < t; ++s) {
			g.nonmove(s, t);
		}
	}
	// render virtual registers' interference
	for (reg t = 0; t < regs; ++t) {
		for (reg s = 0; s < t; ++s) {
			const auto overlap = live.lazy(idx(s, 0)) & live.lazy(idx(t, 0));
			if (bits::reduce_or(overlap).eval(max_index)) {
				g.nonmove(s, t);
			}
		}
	}

	std::cout << "graph:\n" << g << std::endl;

	return g;
}

void strip(graph &interf, reg phys_regs, stack<reg> *stk, bool *stripped)
{
	bool _stripped = false;
	while (true) {
		// find any node of insignificant degree
		// default with an impossible value
		reg chosen_reg = static_cast<reg>(interf.order());
		// else find a node to spill
		// here the heuristic selects the max degree
		// it could be more sophisticated
		size_t max_deg = 0;
		for (reg r = 0; r < interf.order(); ++r) {
			if (!interf.is_nonmove(r))
				continue;
			const auto deg = interf.deg(r);
			if (deg < phys_regs) {
				chosen_reg = r;
				break;
			} else if (deg > max_deg) {
				chosen_reg = r;
				max_deg = deg;
			}
		}
		if (chosen_reg == interf.order()) break;
		stk->push(chosen_reg);
		interf.remove(chosen_reg);
		_stripped = *stripped = true;
	}

	std::cout << "strip:\n";
	if (_stripped)
		std::cout << interf << *stk;
}

void remap(code_t &code, std::vector<reg> const &mapping)
{
	for (auto &ins : code) {
		switch (ins.opcode) {
		case instr::add:
		case instr::umul:
			ins.rs2 = mapping[ins.rs2];
			// fallthrough
		case instr::copy:
		case instr::load:
		case instr::store:
		case instr::clobber:
			ins.rs1 = mapping[ins.rs1];
			// fallthrough
		case instr::def:
		case instr::req:
		case instr::imm:
		case instr::load_local:
		case instr::store_local:
			ins.rd  = mapping[ins.rd ];
			break;
		}
	}
}

void merge(graph &interf, code_t &code, stack<reg> *stk, bool *merged)
{
	bool _merged = false;
	const auto order = interf.order();
	std::vector<reg> mapping(order);
	for (reg r = 0; r < mapping.size(); ++r) {
		mapping[r] = r;
	}
	for (reg s = 0; s < order; ++s) {
		if (!interf.is_live(s))
			continue;
		for (reg t = s + 1; t < order; ++t) {
			// anything >= s is unmapped
			// anything <  s may be mapped
			if (!interf.is_live(t) || !interf.is_move(s, t))
				continue;
			const auto neighbr_nonmove = interf.nonmove_.lazy(interf.idx(s, 0)) | interf.nonmove_.lazy(interf.idx(t, 0));
			const auto neighbr_move    = interf.   move_.lazy(interf.idx(s, 0)) | interf.   move_.lazy(interf.idx(t, 0));
			const auto neighbr         = (neighbr_nonmove | neighbr_move) & ~interf.removed.lazy();
			const auto deg = bits::count(neighbr);

			// Briggs test
			if (deg.eval(order) - 2 /* s and t */ < code.phys_regs) {
				// keep going
			} else {
				// George test
				for (reg unmapped_r = 0; unmapped_r < order; ++unmapped_r) {
					reg r = mapping[unmapped_r];
					if (!interf.is_live(r) || interf.is_free(r, t))
						continue;
					if (interf.is_free(r, s)
					 && interf.deg(r) >= code.phys_regs)
						goto next_iter;
				}
			}
			// small reg takes over large reg
			// with the mapping, t is GONE gone
			mapping[t] = s;
			interf.remove(t);
			neighbr_nonmove.eval(order, interf.nonmove_.slice(interf.idx(s, 0)));
			neighbr_move   .eval(order, interf.   move_.slice(interf.idx(s, 0)));
			interf.move_.clear(interf.idx(s, s));
			_merged = *merged = true;
		next_iter:
			;
		}
	}
	remap(code, mapping);

	std::cout << "merge:\n";
	if (_merged)
		std::cout << interf << *stk << mapping << code;
}

void freeze(graph &interf, reg phys_regs, bool *froze)
{
	bool _froze = false;
	for (reg r = 0; r < interf.order(); ++r) {
		if (!interf.is_live(r))
			continue;
		const auto deg = interf.deg(r);
		if (deg < phys_regs) {
			interf.freeze(r);
			_froze = *froze = true;
			break;
		}
	}

	std::cout << "freeze:\n";
	if (_froze)
		std::cout << interf;
}

void evict(graph &interf, reg phys_regs, stack<reg> *stk, bool *acted)
{
	bool _evicted = false;
	size_t max_deg = 0;
	reg best = interf.order();
	for (reg r = 0; r < interf.order(); ++r) {
		if (!interf.is_live(r))
			continue;
		const auto deg = interf.deg(r);
		if (deg > max_deg) {
			max_deg = deg;
			best = r;
		}
	}
	if (best != interf.order()) {
		stk->push(best);
		interf.remove(best);
		_evicted = *acted = true;
	}

	std::cout << "evict:\n";
	if (_evicted)
		std::cout << interf << *stk;
}

std::vector<reg> select(graph const &interf, stack<reg> stk,
		reg phys_regs, bool *spilled, bits::managed &bound)
{
	auto mapping = std::vector<reg>(interf.order());
	for (reg phys = 0; phys < phys_regs; ++phys) {
		mapping[phys] = phys;
		bound.set(phys);
	}
	while (!stk.empty()) {
		const auto s = stk.pop();
		// represents a mask of free registers relative to neighbors
		reg free = bit(phys_regs) - 1;
		assert(sizeof(free) * CHAR_BIT >= phys_regs);
		if (bound[s])
			continue;
		for (reg t = 0; t < interf.order(); ++t) {
			// clears a bit corresponding to a register that is already mapped
			int msk = !interf.is_free(s, t);
			free &= ~(reg(msk & bound[t]) << mapping[t]);
		}
		if (free) {
			mapping[s] = some_bit_index(free);
			bound.set(s);
			assert(mapping[s] < CHAR_BIT * sizeof(reg));
		} else {
			*spilled = true;
			mapping[s] = s;
		}
	}

	std::cout << "select:\n" << mapping;

	return mapping;
}

code_t rewrite(code_t const &code, bits::managed const &bound)
{
	code_t next_code = { {}, code.phys_regs, code.virt_regs };
	enum usage { use = 0, def = 1 };
	const auto ref = [&] (reg r, usage u) {
		if (r >= code.phys_regs && !bound[r]) {
			const reg next_r = next_code.phys_regs + next_code.virt_regs++;
			static_assert(instr::load_local + 0 == instr:: load_local);
			static_assert(instr::load_local + 1 == instr::store_local);
			const auto opcode = static_cast<instr::opcode_t>(instr::load_local + static_cast<int>(u));
			next_code.emplace_back(opcode, next_r, r);
			r = next_r;
		}
		return r;
	};

	for (auto ins : code) {
		// force patchme to be valid during the iteration
		next_code.reserve(next_code.size() + 4);
		switch (ins.opcode) {
			instr *patchme;

		case instr::copy:
			if (ins.rd != ins.rs1) {
				ins.rs1 = ref(ins.rs1, use);
				patchme = &next_code.emplace_back(ins);
				patchme->rd = ref(ins.rd, def);
			}
			break;

		case instr::add:
		case instr::umul:
			ins.rs2 = ref(ins.rs2, use);
			// fallthrough
		case instr::load:
			ins.rs1 = ref(ins.rs1, use);
			// fallthrough
		case instr::def:
		case instr::imm:
		case instr::load_local:
			patchme = &next_code.emplace_back(ins);
			patchme->rd = ref(ins.rd, def);
			break;

		case instr::store:
			ins.rs1 = ref(ins.rs1, use);
			// fallthrough
		case instr::store_local:
		case instr::req:
			ins.rd  = ref(ins.rd , use);
			next_code.emplace_back(ins);
			break;

		case instr::clobber:
			// a spilled variable isn't going to clobber anyone
			if (bound[ins.rs1])
				next_code.emplace_back(ins);
			break;
		}
	}
	return next_code;
}

std::vector<reg> gcolor(code_t &code)
{
	code = canonical(code);
	std::vector<reg> offsets(code.virt_regs);
	while (true) {
		auto interf = gen_graph(code);
		stack<reg> stk;
		bool acted;
		do {
			do {
				do {
					acted = false;
					strip(interf, code.phys_regs, &stk, &acted);
					merge(interf, code, &stk, &acted);
				} while (acted);
				freeze(interf, code.phys_regs, &acted);
			} while (acted);
			evict(interf, code.phys_regs, &stk, &acted);
		} while (acted);
		bits::managed bound(code.regs());
		const auto color = select(interf, std::move(stk), code.phys_regs, &acted, bound);
		remap(code, color);
		// TODO: apply a spill-less graph coloring to spilled nodes to debloat the graph
		// strip->merge(no degree tests) cycle ->select
		code = rewrite(code, bound);
		if (!acted) {
			return color;
		}
	}
}

int main()
{

	enum { a = 0xa, b, c, d, e, f };
#if 1
	code_t code{
		{
			{ instr::def , 3 },
			{ instr::copy, e, 3 },
			{ instr::def , 0 },
			{ instr::def , 1 },
			{ instr::copy, c, 0 },
			{ instr::copy, b, 1 },
			{ instr::load, 9, b },
			{ instr::add , a, c, 9 },
			{ instr::add , 8, 9, a },
			{ instr::load, 7, b },
			{ instr::load, d, b },
			{ instr::load, 4, 8 },
			{ instr::umul, 5, 8, 4 },
			{ instr::copy, 6, 5 },
			{ instr::add , c, d, 6 },
			{ instr::copy, b, 4 },
			{ instr::copy, 0, 6 },
			{ instr::copy, 1, c },
			{ instr::copy, 2, b },
			{ instr::req , 0 },
			{ instr::req , 1 },
			{ instr::copy, 3, e },
			{ instr::req , 3 },
		},
		4,
		11,
	};
#endif
#if 0
	code_t code = {
		{
			{ instr::def , 0 },
			{ instr::def , 1 },
			{ instr::copy, 3, 0 },
			{ instr::add , 2, 0, 1 },
			{ instr::add , 4, 3, 1 },
			{ instr::add , 5, 2, 4 },
			{ instr::req , 5 },
		},
		3,
		6,
	};
#endif

	std::cout << "original:\n" << code;
	auto m = gcolor(code);
	std::cout << "colored:\n" << code;
	std::cout << "coloring:\n" << m;
	return 0;
}


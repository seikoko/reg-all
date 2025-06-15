#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <deque>
#include <climits>
#include <vector>
#include <cstdint>

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

reg bits(reg r)
{
	return bit(r) - 1;
}

template <typename T, size_t bits = sizeof(T) * CHAR_BIT>
struct bitset
{
	using unit = std::uint32_t;
	enum { unit_bits = sizeof(unit) * CHAR_BIT };
	static_assert(bits <= unit_bits, "large bitfield not supported");

	bitset(unit elems);
	T operator[](unit i) const;
	void set(unit i, T value);
	void acc(unit i, T value);
	void del(unit i, T value);
	unit filter_count(T mask) const;
	bool zero() const;

	std::vector<unit> data;
	unit elems;
};

template <typename T, size_t bits>
bitset<T, bits>::bitset(unit elems)
	: elems(elems)
{
	const auto units = (elems * bits + unit_bits - 1) / unit_bits;
	data.assign(units, 0);
}

// TODO: easily vectorizable
template <typename T, size_t bits>
bitset<T, bits> operator&(bitset<T, bits> const &lhs, bitset<T, bits> const &rhs)
{
	assert(lhs.elems == rhs.elems);
	bitset<T, bits> res(lhs.elems);
	for (size_t i = 0; i < lhs.data.size(); ++i) {
		res.data[i] = lhs.data[i] & rhs.data[i];
	}
	return res;
}

// TODO: easily vectorizable
template <typename T, size_t bits>
bitset<T, bits> operator|(bitset<T, bits> const &lhs, bitset<T, bits> const &rhs)
{
	assert(lhs.elems == rhs.elems);
	bitset<T, bits> res(lhs.elems);
	for (size_t i = 0; i < lhs.data.size(); ++i) {
		res.data[i] = lhs.data[i] | rhs.data[i];
	}
	return res;
}

template <typename T, size_t bits>
T bitset<T, bits>::operator[](std::uint32_t i) const
{
	assert(i < elems);
	unit res = 0;
	for (size_t b = 0; b < bits; ++b) {
		size_t bit_idx = i + b * elems;
		res |= (data[bit_idx / unit_bits] >> (bit_idx % unit_bits) & 1) << b;
	}
	return static_cast<T>(res);
}

template <typename T, size_t bits>
void bitset<T, bits>::set(unit i, T _value)
{
	assert(i < elems);
	unit value = _value;
	for (size_t b = 0; b < bits; ++b) {
		size_t bit_idx = i + b * elems;
		data[bit_idx / unit_bits] &= ~(unit{1} << (bit_idx % unit_bits));
		data[bit_idx / unit_bits] |= (value & 1) << (bit_idx % unit_bits);
		value >>= 1;
	}
}

template <typename T, size_t bits>
void bitset<T, bits>::acc(unit i, T _value)
{
	assert(i < elems);
	unit value = _value;
	for (size_t b = 0; b < bits; ++b) {
		size_t bit_idx = i + b * elems;
		data[bit_idx / unit_bits] |= (value & 1) << (bit_idx % unit_bits);
		value >>= 1;
	}
}

template <typename T, size_t bits>
void bitset<T, bits>::del(unit i, T _value)
{
	assert(i < elems);
	unit value = _value;
	for (size_t b = 0; b < bits; ++b) {
		size_t bit_idx = i + b * elems;
		data[bit_idx / unit_bits] &= ~((value & 1) << (bit_idx % unit_bits));
		value >>= 1;
	}
}

template <typename T, size_t bits>
typename bitset<T, bits>::unit bitset<T, bits>::filter_count(T mask) const
{
	unit count = 0;
	for (size_t i = 0; i < data.size(); ++i) {
		unit filter = 0;
		for (reg b = 0; b < bits; ++b) {
			if (bit(b) & mask)
				filter |= data[i + b * elems];
		}
		count += __builtin_popcount(filter);
	}
	return count;
}

template <typename T, size_t bits>
bool bitset<T, bits>::zero() const
{
	for (size_t i = 0; i < data.size(); ++i) {
		if (data[i])
			return false;
	}
	return true;
}

struct graph
{
	enum interference { I_FREE = 0, I_MOVE = 1, I_NONMOVE = 2 /* or 3 */ };
	static bool is_free(interference i) { return i == I_FREE; }
	static bool is_move(interference i) { return i == I_MOVE; }
	static bool is_nonmove(interference i) { return i >= I_NONMOVE; }

	std::vector<bitset<interference, 2>> madj;

	enum state { S_VISIBLE = 0, S_FROZEN = 1, S_HIDDEN = 2 };
	bitset<state, 2> removed;

	graph(reg count);

	void link(reg s, reg t, interference i);
	interference linked(reg s, reg t) const;
	void remove(reg s);
	bool has(reg s) const;
	reg order() const;
	reg deg(reg s) const;
	bool nonmove(reg s) const;
	reg move_related(reg s) const;
};

graph::graph(reg count)
	: madj(count, count), removed(count)
{
}

reg graph::order() const
{
	return removed.elems;
}

void graph::link(reg s, reg t, interference i)
{
	// assert(i != NONMOVE || linked(s, t) != I_NONMOVE);
	madj[s].acc(t, i);
	madj[t].acc(s, i);
}

graph::interference graph::linked(reg s, reg t) const
{
	return madj[s][t];
}

void graph::remove(reg s)
{
	// assert(removed[s] != S_HIDDEN);
	removed.set(s, S_HIDDEN);
}

bool graph::has(reg s) const
{
	return removed[s] == S_VISIBLE;
}

reg graph::deg(reg s) const
{
	reg d = 0;
	for (reg t = 0; t < order(); ++t) {
		if (has(t) && linked(s, t))
			++d;
	}
	return d;
}

reg graph::move_related(reg s) const
{
	reg m = 0;
	for (reg t = 0; t < order(); ++t) {
		if (has(t) && is_move(linked(s, t)))
			++m;
	}
	return m;
}

bool graph::nonmove(reg s) const
{
	return removed[s] == S_FROZEN || (removed[s] == S_VISIBLE && !move_related(s));
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

template <typename T, size_t bits>
std::ostream &operator<<(std::ostream &os, bitset<T, bits> const &b)
{
	os << "| ";
	for (typename bitset<T, bits>::unit i = 0; i < b.elems; ++i) {
		if (b[i]) {
			os << i;
			if (bits > 1)
				os << "." << b[i];
			os << " ";
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
		if (g.removed[r] & graph::S_HIDDEN)
			continue;
		os << " F"[g.removed[r]];
		os << r << "(d=" << g.deg(r) << ",m=" << g.move_related(r) << "): ";
		bool first = true;
		for (reg t = 0; t < g.order(); ++t) {
			if (!g.linked(r, t) || !g.has(t))
				continue;
			if (first) {
				first = false;
			} else {
				os << " ";
			}
			os << " M ."[g.linked(r, t)] << t;
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
	graph g(code.regs());
	// [definition, last use)
	std::vector<bitset<bool, 1>> live(code.regs(), code.size() + 1);
	const auto use = [&] (reg r, reg index) { live[r].acc(index, true); };
	const auto def = [&] (reg r, reg index) { live[r].del(index, true); };
	const auto clobber = [&] (reg virt, reg phys) { g.link(virt, phys, graph::I_NONMOVE); };

	static int _iter = 0;
	std::cout << "iter #" << _iter++ << ":\n";

	reg i = code.size() - 1;
	do {
		for (reg r = 0; r < code.regs(); ++r) {
			live[r].set(i, live[r][i+1]);
		}

		const auto ins = code[i];
		switch (ins.opcode) {
		case instr::copy:
			def(ins.rd , i);
			use(ins.rs1, i);
			g.link(ins.rd, ins.rs1, graph::I_MOVE);
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
			def(ins.rd, i);
			use(ins.rs1, i);
			break;
		case instr::store_local:
			use(ins.rd, i);
			break;
		case instr::def:
			def(ins.rd, 0);
			break;
		case instr::req:
			use(ins.rd, code.size()-1);
			break;
		case instr::clobber:
			clobber(ins.rd, ins.rs1);
			break;
		}
	} while (i--);

	for (reg i = 0; i < code.size(); ++i) {
		bitset<bool, 1> _live(code.regs());
		for (reg r = 0; r < code.regs(); ++r) {
			if (live[r][i])
				_live.acc(r, true);
		}
		std::cout << i << ": " << _live << code[i];
	}

	// render physical registers
	for (reg t = 1; t < code.phys_regs; ++t) {
		for (reg s = 0; s < t; ++s) {
			g.link(s, t, graph::I_NONMOVE);
		}
	}
	// render virtual registers' interference
	for (reg t = 0; t < code.regs(); ++t) {
		for (reg s = 0; s < t; ++s) {
			const auto overlap = live[s] & live[t];
			if (!overlap.zero()) {
				g.link(s, t, graph::I_NONMOVE);
			}
		}
	}

	std::cout << "graph:\n" << g << std::endl;

	return g;
}

void strip(graph &interference, reg phys_regs, stack<reg> *stk, bool *stripped)
{
	bool _stripped = false;
	while (true) {
		// find any node of insignificant degree
		// default with an impossible value
		reg chosen_reg = static_cast<reg>(interference.order());
		// else find a node to spill
		// here the heuristic selects the max degree
		// it could be more sophisticated
		reg max_deg = 0;
		for (reg r = 0; r < interference.order(); ++r) {
			if (!interference.nonmove(r))
				continue;
			const auto deg = interference.deg(r);
			if (deg < phys_regs) {
				chosen_reg = r;
				break;
			} else if (deg > max_deg) {
				chosen_reg = r;
				max_deg = deg;
			}
		}
		if (chosen_reg == interference.order()) break;
		stk->push(chosen_reg);
		// TODO: potential optimization : just set interference.deg(r) = virt_regs+1 so removed nodes are never selected
		interference.remove(chosen_reg);
		_stripped = *stripped = true;
	}

	std::cout << "strip:\n";
	if (_stripped)
		std::cout << interference << *stk;
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

void merge(graph &interference, code_t &code, stack<reg> *stk, bool *merged)
{
	bool _merged = false;
	std::vector<reg> mapping(interference.order());
	for (reg r = 0; r < mapping.size(); ++r) {
		mapping[r] = r;
	}
	for (reg s = 0; s < interference.order(); ++s) {
		if (!interference.has(s))
			continue;
		for (reg t = s + 1; t < interference.order(); ++t) {
			// anything >= s is unmapped
			// anything <  s may be mapped
			if (!interference.has(t) || !graph::is_move(interference.linked(s, t)))
				continue;
			auto neighbours = interference.madj[s] | interference.madj[t];
			for (reg unmapped_r = 0; unmapped_r < interference.order(); ++unmapped_r) {
				reg r = mapping[unmapped_r];
				if (!interference.has(r))
					neighbours.set(r, graph::I_FREE);
			}
			neighbours.set(s, graph::I_FREE);
			neighbours.set(t, graph::I_FREE);
			// I_MOVE | I_NONMOVE
			const auto deg = neighbours.filter_count(graph::interference(3));
			// Briggs test
			if (deg < code.phys_regs) {
				// keep going
			} else {
				// George test
				for (reg unmapped_r = 0; unmapped_r < interference.order(); ++unmapped_r) {
					reg r = mapping[unmapped_r];
					assert(!interference.has(r) || (interference.linked(r, t) == interference.linked(t, r)));
					if (!interference.has(r) || graph::is_free(interference.linked(r, t)))
						continue;
					if (graph::is_free(interference.linked(r, s))
					 && interference.deg(r) >= code.phys_regs)
						goto next_iter;
				}
			}
			// small reg takes over large reg
			mapping[t] = s;
			interference.remove(t);
			stk->push(t);
			for (reg unmapped_r = 0; unmapped_r < interference.order(); ++unmapped_r) {
				reg r = mapping[unmapped_r];
				if (neighbours[r])
					interference.link(r, s, neighbours[r]);
			}
			interference.madj[s] = neighbours;
			_merged = *merged = true;
		next_iter:
			;
		}
	}
	remap(code, mapping);

	std::cout << "merge:\n";
	if (_merged)
		std::cout << interference << *stk << mapping << code;
}

void freeze(graph &interference, reg phys_regs, bool *froze)
{
	bool _froze = false;
	for (reg r = 0; r < interference.order(); ++r) {
		if (!interference.has(r))
			continue;
		const auto deg = interference.deg(r);
		if (deg < phys_regs) {
			interference.removed.acc(r, graph::S_FROZEN);
			_froze = *froze = true;
			break;
		}
	}

	std::cout << "freeze:\n";
	if (_froze)
		std::cout << interference;
}

void evict(graph &interference, reg phys_regs, stack<reg> *stk, bool *acted)
{
	bool _evicted = false;
	reg max_deg = 0;
	reg best = interference.order();
	for (reg r = 0; r < interference.order(); ++r) {
		if (interference.removed[r] == graph::S_HIDDEN)
			continue;
		const auto deg = interference.deg(r);
		if (deg > max_deg) {
			max_deg = deg;
			best = r;
		}
	}
	if (best != interference.order()) {
		stk->push(best);
		interference.remove(best);
		_evicted = *acted = true;
	}

	std::cout << "evict:\n";
	if (_evicted)
		std::cout << interference << *stk;
}

std::vector<reg> select(graph const &interference, stack<reg> stk,
		reg phys_regs, bool *spilled, bitset<bool, 1> &bound)
{
	auto mapping = std::vector<reg>(interference.order());
	for (reg phys = 0; phys < phys_regs; ++phys) {
		mapping[phys] = phys;
		bound.acc(phys, true);
	}
	while (!stk.empty()) {
		const auto s = stk.pop();
		// represents a mask of free registers relative to neighbors
		reg free = bits(phys_regs);
		assert(sizeof(free) * CHAR_BIT >= phys_regs);
		if (bound[s])
			continue;
		for (reg t = 0; t < interference.order(); ++t) {
			// clears a bit corresponding to a register that is already mapped
			int msk = interference.madj[s][t];
			msk |= msk >> 1;
			free &= ~(reg(msk & bound[t]) << mapping[t]);
		}
		if (free) {
			mapping[s] = some_bit_index(free);
			bound.acc(s, true);
			assert(mapping[s] < CHAR_BIT * sizeof(reg));
		} else {
			*spilled = true;
			mapping[s] = s;
		}
	}

	std::cout << "select:\n" << mapping;

	return mapping;
}

code_t rewrite(code_t const &code, bitset<bool, 1> const &bound)
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
		auto interference = gen_graph(code);
		stack<reg> stk;
		bool acted;
		do {
			do {
				do {
					acted = false;
					strip(interference, code.phys_regs, &stk, &acted);
					merge(interference, code, &stk, &acted);
				} while (acted);
				freeze(interference, code.phys_regs, &acted);
			} while (acted);
			evict(interference, code.phys_regs, &stk, &acted);
		} while (acted);
		assert(stk.size() >= code.virt_regs);
		bitset<bool, 1> bound(code.regs());
		const auto color = select(interference, std::move(stk), code.phys_regs, &acted, bound);
		remap(code, color);
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
			// { instr::def , 3 },
			// { instr::copy, e, 3 },
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
			// { instr::copy, 1, c },
			// { instr::copy, 2, b },
			{ instr::req , 0 },
			// { instr::req , 1 },
			// { instr::copy, 3, e },
			// { instr::req , 3 },
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


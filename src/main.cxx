#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <deque>
#include <climits>
#include <vector>
#include <cstdint>

using reg = std::uint32_t;

struct instr
{
	enum opcode_t : reg {
		def,   // rd <- ?
		req,   // rd -> ?
		copy,  // rd <- rs1
		load,  // rd <- [rs1]
		store, // rd -> [rs1]
		load_local,  // rd <- local_#rs1
		store_local, // rd -> local_#rs1
		add,   // rd <- rs1 + rs2
		imm,   // rd <- #rs1
	} opcode;
	reg rd ;
	reg rs1;
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
	static_assert(bits != 0 && (bits & (bits-1)) == 0, "non power of 2 bitcount not supported");

	bitset(unit elems);
	T operator[](unit i) const;
	void set(unit i, T value);

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

template <typename T, size_t bits>
T bitset<T, bits>::operator[](std::uint32_t i) const
{
	assert(i < elems);
	return static_cast<T>((data[i * bits / unit_bits] >> (i * bits % unit_bits)) & ::bits(bits));
}

template <typename T, size_t bits>
void bitset<T, bits>::set(unit i, T value)
{
	assert(i < elems);
	unit msk = ::bits(bits) << (i * bits % unit_bits);
	data[i * bits / unit_bits] &= ~msk;
	data[i * bits / unit_bits] |=  msk & static_cast<unit>(value << (i * bits % unit_bits));
}

struct graph
{
	std::vector<std::vector<reg>> ladj;
	std::vector<reg> degree;
	bitset<bool, 1> removed;

	graph(reg count);

	void link(reg s, reg t);
	void remove(reg s);
	bool has(reg s) const;
	size_t order() const;
};

graph::graph(reg count)
	: ladj(count), degree(count), removed(count)
{
}

size_t graph::order() const
{
	assert(degree.size() == ladj.size());
	assert(removed.elems == ladj.size());
	return ladj.size();
}

void graph::link(reg s, reg t)
{
	if (s > t) std::swap(s, t);
	ladj[s].push_back(t);
	ladj[t].push_back(s);
	degree[s]++;
	degree[t]++;
}

void graph::remove(reg s)
{
	removed.set(s, true);
	for (const auto t : ladj[s]) {
		degree[t]--;
	}
}

bool graph::has(reg s) const
{
	return !removed[s];
}

template <typename T>
std::ostream &operator<<(std::ostream &os, std::vector<T> const &v)
{
	os << "{\n";
	for (size_t i = 0; i < v.size(); ++i) {
		os << i << ": " << v[i] << '\n';
	}
	os << "}\n";
	return os;
}

std::ostream &operator<<(std::ostream &os, instr const &ins)
{
	switch (ins.opcode) {
	case instr::def:
		return os << "def   %" << ins.rd;
	case instr::req:
		return os << "req   %" << ins.rd;
	case instr::copy:
		return os << "copy  %" << ins.rd << ", %" << ins.rs1;
	case instr::load:
		return os << "load  %" << ins.rd << ", [%" << ins.rs1 << "]";
	case instr::store:
		return os << "store %" << ins.rd << ", [%" << ins.rs1 << "]";
	case instr::add:
		return os << "add   %" << ins.rd << ", %" << ins.rs1 << ", %" << ins.rs2;
	case instr::imm:
		return os << "imm   %" << ins.rd << ", #" << ins.rs1;
	case instr::load_local:
		return os << "load  %" << ins.rd << ", local_" << ins.rs1;
	case instr::store_local:
		return os << "store %" << ins.rd << ", local_" << ins.rs1;
	default:
		assert(false);
	}
}

std::ostream &operator<<(std::ostream &os, graph const &g)
{
	for (reg r = 0; r < g.order(); ++r) {
		os << r << "(" << g.degree[r] << "): ";
		bool first = true;
		if (g.has(r)) for (auto t : g.ladj[r]) {
			if (!g.has(t))
				continue;
			if (first) {
				first = false;
			} else {
				os << ", ";
			}
			os << t;
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

graph gen_graph(code_t const &code)
{
	// first registers correspond to the physical registers and all interfere with each other
	// currently instructions must not refer to physical registers directly
	// TODO: allow instructions to refer to phys registers
	// enabling `mul v0, v1, v2;`
	// to transform into `clobber v0, edx; copy eax, v1; mul eax, eax, v2; copy v0, eax`
	graph g(code.regs());
	// [definition, last use)
	std::vector<std::pair<reg, reg>> live;
	for (reg r = 0; r < code.virt_regs; ++r) {
		live.emplace_back(reg(-1), reg(0));
	}
	const auto idx = [&] (reg i) -> auto&   { return assert(i >= code.phys_regs), live[i - code.phys_regs]; };
	const auto use = [&] (reg r, reg index) { if (idx(r).second < index) idx(r).second = index; };
	const auto def = [&] (reg r, reg index) { if (idx(r).first  > index) idx(r).first  = index; };
	const auto clobber = [&] (reg virt, reg phys) { g.link(virt, phys); };

	for (reg i = 0; i < code.size(); ++i) {
		const auto ins = code[i];
		switch (ins.opcode) {
		case instr::add:
			use(ins.rs2, i);
			clobber(ins.rd, 0);
			// fallthrough
		case instr::copy:
		case instr::load:
			use(ins.rs1, i);
			// fallthrough
		case instr::def:
		case instr::imm:
		case instr::load_local:
			use(ins.rd, i);
			def(ins.rd, i);
			break;
		case instr::store:
			use(ins.rs1, i);
			// fallthrough
		case instr::store_local:
			use(ins.rd, i);
			break;
		case instr::req:
			use(ins.rd, i);
			break;
		}

	}
	for (reg r = 0; r < code.virt_regs; ++r) {
		assert(live[r].second != 0 || live[r].first == reg(-1));
	}

	// render physical registers
	for (reg t = 1; t < code.phys_regs; ++t) {
		for (reg s = 0; s < t; ++s) {
			g.link(s, t);
		}
	}
	// render virtual registers' interference
	for (reg t = code.phys_regs + 1; t < code.regs(); ++t) {
		for (reg s = code.phys_regs; s < t; ++s) {
			bool overlap = idx(s).first < idx(t).second && idx(t).first < idx(s).second;
			if (overlap) {
				g.link(s, t);
			}
		}
	}
	return g;
}

stack<reg> strip(graph &interference, reg phys_regs)
{
	stack<reg> stk;
	while (true) {
		// find any node of insignificant degree
		// default with an impossible value
		reg chosen_reg = static_cast<reg>(interference.order());
		// else find a node to spill
		// here the heuristic selects the max degree
		// it could be more sophisticated
		reg max_deg = 0;
		for (reg r = phys_regs; r < interference.order(); ++r) {
			if (!interference.has(r))
				continue;
			const auto deg = interference.degree[r];
			if (deg < phys_regs) {
				chosen_reg = r;
				break;
			} else if (deg > max_deg) {
				chosen_reg = r;
				max_deg = deg;
			}
		}
		if (chosen_reg == interference.order()) break;
		stk.push(chosen_reg);
		// TODO: potential optimization : just set interference.degree[r] = virt_regs+1 so removed nodes are never selected
		interference.remove(chosen_reg);
	}
	return stk;
}

std::vector<reg> select(graph const &interference, stack<reg> stk, reg phys_regs, bool *spilled, bitset<bool, 1> &bound)
{
	auto mapping = std::vector<reg>(interference.order());
	for (reg phys = 0; phys < phys_regs; ++phys) {
		mapping[phys] = phys;
		bound.set(phys, true);
	}
	while (!stk.empty()) {
		const auto s = stk.pop();
		// represents a mask of free neighboring registers
		reg free = bits(phys_regs);
		assert(sizeof(free) * CHAR_BIT >= phys_regs);
		for (const auto t : interference.ladj[s]) {
			// clears a bit corresponding to a register that is already mapped
			free &= ~(reg(bound[t]) << mapping[t]);
		}
		bound.set(s, !!free);
		if (free) {
			mapping[s] = some_bit_index(free);
			assert(mapping[s] < CHAR_BIT * sizeof(reg));
		} else {
			*spilled = true;
		}
	}
	return mapping;
}

code_t rewrite(code_t const &code, bitset<bool, 1> const &bound)
{
	code_t next_code = { {}, code.phys_regs, code.virt_regs };
	enum usage { use, def };
	const auto ref = [&] (reg r, usage u) {
		assert(r >= code.phys_regs);
		if (!bound[r]) {
			const reg next_r = next_code.phys_regs + next_code.virt_regs++;
			static_assert(instr::load_local + 1 == instr::store_local);
			const auto opcode = static_cast<instr::opcode_t>(instr::load_local + (u == def));
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
		case instr::def:
		case instr::imm:
			patchme = &next_code.emplace_back(ins);
			patchme->rd = ref(ins.rd, def);
			break;
		case instr::req:
			ins.rd = ref(ins.rd, use);
			next_code.emplace_back(ins);
			break;
		case instr::copy:
		case instr::load:
			ins.rs1 = ref(ins.rs1, use);
			patchme = &next_code.emplace_back(ins);
			patchme->rd = ref(ins.rd, def);
			break;
		case instr::add:
			ins.rs1 = ref(ins.rs1, use);
			ins.rs2 = ref(ins.rs2, use);
			patchme = &next_code.emplace_back(ins);
			patchme->rd = ref(ins.rd, def);
			break;
		case instr::store:
			ins.rs1 = ref(ins.rs1, use);
			ins.rd  = ref(ins.rd , use);
			next_code.emplace_back(ins);
			break;
		case instr::load_local:
			// maybe could delete the variable if this is executed
			patchme = &next_code.emplace_back(ins);
			patchme->rd = ref(ins.rd, def);
			break;
		case instr::store_local:
			ins.rd = ref(ins.rd, use);
			next_code.emplace_back(ins);
			break;
		}
	}
	return next_code;
}

void remap(code_t &code, std::vector<reg> const &mapping)
{
	for (auto &ins : code) {
		switch (ins.opcode) {
		case instr::add:
			ins.rs2 = mapping[ins.rs2];
			// fallthrough
		case instr::copy:
		case instr::load:
		case instr::store:
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

std::vector<reg> gcolor(code_t &code)
{
	std::vector<reg> offsets(code.virt_regs);
	for (reg r = 0; r < offsets.size(); ++r) {
		offsets[r] = r + code.phys_regs;
	}
	remap(code, offsets);
	while (true) {
		auto interference = gen_graph(code);
		stack<reg> stk = strip(interference, code.phys_regs);
		bool spilled = false;
		bitset<bool, 1> bound(code.regs());
		const auto mapping = select(interference, std::move(stk), code.phys_regs, &spilled, bound);
		code = rewrite(code, bound);
		if (!spilled) {
			remap(code, mapping);
			return mapping;
		}
	}
}

int main()
{
	code_t code{
		{
			{ instr::def , 8 },
			{ instr::def , 7 },
			{ instr::load, 5, 7 },
			{ instr::add , 6, 8, 5 },
			{ instr::add , 4, 5, 6 },
			{ instr::load, 3, 7 },
			{ instr::load, 9, 7 },
			{ instr::load, 0, 4 },
			{ instr::add , 1, 3, 0 },
			{ instr::copy, 2, 1 },
			{ instr::add , 8, 9, 2 },
			{ instr::copy, 7, 0 },
			{ instr::req , 2 },
			{ instr::req , 8 },
			{ instr::req , 7 },
		},
		4,
		10,
	};
	std::cout << code;
	auto m = gcolor(code);
	std::cout << code;
	std::cout << m;
	return 0;
}


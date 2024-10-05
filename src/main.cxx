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
	return reg{1} << (r % (CHAR_BIT * sizeof r));
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
};

graph::graph(reg count)
	: ladj(count), degree(count), removed(count)
{
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

graph gen_graph(reg count, std::vector<instr> const &v)
{
	graph g(count);
	// [definition, last use)
	std::vector<std::pair<reg, reg>> live;
	for (reg r = 0; r < count; ++r) {
		live.emplace_back(reg(-1), reg(0));
	}
	const auto use = [&] (reg r, reg index) { if (live[r].second < index) live[r].second = index; };
	const auto def = [&] (reg r, reg index) { if (live[r].first  > index) live[r].first  = index; };
	for (reg i = 0; i < v.size(); ++i) {
		const auto ins = v[i];
		switch (ins.opcode) {
		case instr::add:
			use(ins.rs2, i);
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
	for (reg r = 0; r < count; ++r) {
		assert(live[r].second != 0 || live[r].first == reg(-1));
	}
	for (reg t = 1; t < count; ++t) {
		for (reg s = 0; s < t; ++s) {
			bool overlap = live[s].first < live[t].second && live[t].first < live[s].second;
			if (overlap) {
				g.link(s, t);
			}
		}
	}
	return g;
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
	for (reg r = 0; r < g.ladj.size(); ++r) {
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

stack<reg> strip(graph &interference, reg phys_regs, reg virt_regs)
{
	stack<reg> stk;
	while (true) {
		// find any node of insignificant degree
		reg chosen_reg = virt_regs + 1;
		for (reg r = 0; r < virt_regs; ++r) {
			if (!interference.has(r))
				continue;
			if (interference.degree[r] < phys_regs) {
				chosen_reg = r;
				break;
			}
		}

		if (chosen_reg == phys_regs + 1) break;
		stk.push(chosen_reg);
		// TODO: potential optimization : just set interference.degree[r] = virt_regs+1 so removed nodes are never selected
		interference.remove(chosen_reg);
	}
	return stk;
}

enum region
{
	virtual_register,
	physical_register,
	potential_spill,
	actual_spill,
};

std::vector<reg> assign(graph const &interference, stack<reg> stk, reg virt_regs, reg phys_regs, bool *spilled, bitset<region, 2> &assignment)
{
	auto mapping = std::vector<reg>(virt_regs);
	while (!stk.empty()) {
		const auto s = stk.pop();
		// represents a mask of neighboring registers in use
		// TODO: should be a bitset
		reg used = 0;
		for (const auto t : interference.ladj[s]) {
			if (assignment[t] == region::physical_register)
				used |= bit(mapping[t]);
		}
		// equivalent of bitwise NOT when you only consider the phys_regs first bits
		const reg free = bits(phys_regs) ^ used;
		if (free) {
			assignment.set(s, region::physical_register);
			// TODO: use some form of heuristic
			mapping[s] = some_bit_index(free);
		} else {
			assignment.set(s, region::actual_spill);
			*spilled = true;
		}
	}
	return mapping;
}

std::vector<instr> rewrite(std::vector<instr> const &code, reg virt_regs, reg *pnext_virt_regs, bitset<region, 2> const &assignment)
{
	auto next_virt_regs = virt_regs;
	std::vector<instr> next_code;

	enum usage { use, def };
	const auto ref = [&] (reg r, usage u) {
		if (assignment[r] == region::actual_spill) {
			const reg next_r = next_virt_regs++;
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
	*pnext_virt_regs = next_virt_regs;
	return next_code;
}

std::vector<reg> gcolor(reg phys_regs, reg virt_regs, std::vector<instr> &code)
{
	while (true) {
		auto interference = gen_graph(virt_regs, code);
		stack<reg> stk = strip(interference, phys_regs, virt_regs);
		for (reg r = 0; r < virt_regs; ++r) {
			if (interference.has(r)) {
				stk.push(r);
			}
		}
		bool spilled = false;
		bitset<region, 2> assignment(virt_regs);
		const auto mapping = assign(interference, std::move(stk), virt_regs, phys_regs, &spilled, assignment);
		code = rewrite(code, virt_regs, &virt_regs, assignment);
		if (!spilled)
			return mapping;
	}
}

int main()
{
	std::vector<instr> v{
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
	};
	std::cout << v;
	auto g = gen_graph(10, v);
	std::cout << g;
	auto m = gcolor(4, 10, v);
	std::cout << v;
	std::cout << m;
	return 0;
}


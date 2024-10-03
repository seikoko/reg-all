#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <deque>
#include <climits>
#include <vector>
#include <cstdint>

using reg = std::uint32_t;

struct color
{
	enum { virt_reg, potential_spill, actual_spill, phys_reg } status;
	reg address;
};

std::ostream &operator<<(std::ostream &os, color c) { return os << "vpsr"[c.status] << std::hex << c.address; }

struct instr
{
	enum {
		def,   // rd <- ?
		req,   // rd -> ?
		copy,  // rd <- rs1
		load,  // rd <- [rs1]
		store, // rd -> [rs1]
		add,   // rd <- rs1 + rs2
		imm,   // rd <- value(rs1)
	} opcode;
	reg rd ;
	reg rs1;
	reg rs2;
};

struct bitset
{
	using unit = std::uint32_t;
	enum { unit_bits = sizeof(unit) * 8 };

	bitset(unit bits);
	bool operator[](unit i) const;
	void set(unit i);
	void unset(unit i);
	void push();

	std::vector<unit> data;
	unit bits;
};

bitset::bitset(unit bits)
	: bits(bits)
{
	const auto units = (bits + unit_bits - 1) / unit_bits;
	data.assign(units, 0);
}

bool bitset::operator[](std::uint32_t i) const
{
	assert(i < bits);
	return (data[i / unit_bits] >> (i % unit_bits)) & 1;
}

void bitset::set(unit i)
{
	assert(i < bits);
	unit msk = unit{1} << (i % unit_bits);
	data[i / unit_bits] &= ~msk;
	data[i / unit_bits] |=  msk;
}

void bitset::unset(unit i)
{
	assert(i < bits);
	unit msk = unit{1} << (i % unit_bits);
	data[i / unit_bits] &= ~msk;
}

void bitset::push()
{
	const auto units = (++bits + unit_bits - 1) / unit_bits;
	if (units > data.size()) {
		data.push_back(0);
	}
}

struct graph
{
	std::vector<std::vector<reg>> ladj;
	std::vector<reg> degree;
	bitset removed;

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
	removed.set(s);
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
	std::vector<std::pair<reg, reg>> live;
	for (reg r = 0; r < count; ++r) {
		live.emplace_back(reg(-1), reg(0));
	}
	for (reg i = 0; i < v.size(); ++i) {
		const auto ins = v[i];
		if (ins.opcode != instr::store) {
			const reg defined = ins.rd;
			if (live[defined].first > i) {
				live[defined].first = i;
			}
		}
		switch (ins.opcode) {
		case instr::add:
			if (live[ins.rs2].second < i) {
				live[ins.rs2].second = i;
			}
			// fallthrough
		case instr::copy:
		case instr::load:
		case instr::store:
			if (live[ins.rs1].second < i) {
				live[ins.rs1].second = i;
			}
			// fallthrough
		case instr::def:
		case instr::imm:
			break;
		case instr::req:
			if (live[ins.rd].second < i) {
				live[ins.rd].second = i;
			}
			break;
		}
	}
	for (reg r = 0; r < count; ++r) {
		assert(live[r].second != 0);
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
		os << "def   %" << ins.rd;
		break;
	case instr::req:
		os << "req   %" << ins.rd;
		break;
	case instr::copy:
		os << "copy  %" << ins.rd << ", %" << ins.rs1;
		break;
	case instr::load:
		os << "load  %" << ins.rd << ", [%" << ins.rs1 << "]";
		break;
	case instr::store:
		os << "store %" << ins.rd << ", [%" << ins.rs1 << "]";
		break;
	case instr::add:
		os << "add   %" << ins.rd << ", %" << ins.rs1 << ", %" << ins.rs2;
		break;
	case instr::imm:
		os << "imm   %" << ins.rd << ", #" << ins.rs1;
		break;
	}
	return os;
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
	return reg{1} << (r % (CHAR_BIT * sizeof r));
}

reg bits(reg r)
{
	return bit(r) - 1;
}

stack<reg> strip(graph &g, reg n_reg, reg v_reg)
{
	stack<reg> stk;
	for (reg best, best_reg;;) {
		// find min
		best = n_reg;
		best_reg = v_reg + 1;
		for (reg r = 0; r < v_reg; ++r) {
			if (!g.has(r))
				continue;
			if (g.degree[r] < best) {
				best = g.degree[r];
				best_reg = r;
			}
		}

		if (best_reg == v_reg + 1) break;
		stk.push(best_reg);
		// TODO: potential optimization : just set g.degree[r] = v_reg+1 so removed nodes are never selected
		g.remove(best_reg);
	}
	return stk;
}

std::vector<color> gcolor(reg n_reg, reg v_reg, std::vector<instr> const &v)
{
	auto g = gen_graph(v_reg, v);
	std::vector<color> mapping(v_reg);
	stack<reg> stk = strip(g, n_reg, v_reg);

	for (reg r = 0; r < v_reg; ++r) {
		if (g.has(r)) {
			mapping[r].status = color::potential_spill;
			mapping[r].address = r;
			g.remove(r);
			stk.push(r);
		}
	}

	while (!stk.empty()) {
		const auto s = stk.pop();
		// represents a mask of neighboring registers in use
		reg used = 0;
		for (const auto t : g.ladj[s]) {
			if (mapping[t].status == color::phys_reg)
				used |= reg{1} << mapping[t].address;
		}
		// equivalent of bitwise NOT when you only consider the n_reg first bits
		const reg free = bits(n_reg) ^ used;
		if (free) {
			mapping[s].status  = color::phys_reg;
			mapping[s].address = some_bit_index(free);
		} else {
			mapping[s].status = color::actual_spill;
		}
	}

	return mapping;
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
	std::cout << m;
	return 0;
}


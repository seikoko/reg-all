#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <deque>
#include <iterator>
#include <vector>
#include <cstdint>

using reg = std::uint32_t;

struct color
{
	enum { virt_reg, potential_spill, actual_spill, phys_reg } status;
	reg address;
};

std::ostream &operator<<(std::ostream &os, color c) { return os << "VPSR"[c.status] << std::hex << c.address; }

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
	void revive(reg s);
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

void graph::revive(reg s)
{
	removed.unset(s);
	for (const auto t : ladj[s]) {
		degree[t]++;
	}
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
	for (reg t = 1; t < count; ++t) {
		for (reg s = 0; s < t; ++s) {
			bool overlap = live[s].first <= live[t].second && live[t].first <= live[s].second;
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
	for (const auto &elem : v) {
		os << elem << '\n';
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
		if (!g.removed[r]) for (auto t : g.ladj[r]) {
			if (g.removed[t])
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

void gcolor(reg count, std::vector<instr> const &v)
{
	const reg n_reg = 6;
	auto g = gen_graph(count, v);
	stack<reg> stk;
	for (reg best, best_reg;;) {
		// find min
		best = n_reg;
		best_reg = count + 1;
		for (reg r = 0; r < count; ++r) {
			if (g.removed[r])
				continue;
			if (g.degree[r] < best) {
				best = g.degree[r];
				best_reg = r;
			}
		}

		if (best_reg == count + 1) break;
		stk.push(best_reg);
		g.remove(best_reg);
	}

	std::cout << g;

	std::vector<color> mapping(count);
	while (!stk.empty()) {
		const auto s = stk.pop();
		reg used = 0;
		for (const auto t : g.ladj[s]) {
			if (g.removed[t])
				continue;
			if (mapping[t].status == color::phys_reg)
				used |= reg{1} << mapping[t].address;
		}
		reg bit_free = lsb(~used);
		mapping[s].status = color::phys_reg;
		mapping[s].address = ctz(bit_free);
		g.revive(s);
	}

	std::cout << mapping;
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
	gcolor(10, v);
	return 0;
}


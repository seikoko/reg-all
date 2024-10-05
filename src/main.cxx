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
		none, // illegal
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
		case instr::none:
			assert(false);
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

stack<reg> strip(graph &interference, reg n_reg, reg v_reg)
{
	stack<reg> stk;
	for (reg best, best_reg;;) {
		// find min
		best = n_reg;
		best_reg = v_reg + 1;
		for (reg r = 0; r < v_reg; ++r) {
			if (!interference.has(r))
				continue;
			if (interference.degree[r] < best) {
				best = interference.degree[r];
				best_reg = r;
			}
		}

		if (best_reg == v_reg + 1) break;
		stk.push(best_reg);
		// TODO: potential optimization : just set interference.degree[r] = v_reg+1 so removed nodes are never selected
		interference.remove(best_reg);
	}
	return stk;
}

std::vector<color> gcolor(reg n_reg, reg v_reg, std::vector<instr> &code)
{
	bool spilled;
	std::vector<color> mapping;
	do {
		spilled = false;
		mapping = std::vector<color>(v_reg);
		auto interference = gen_graph(v_reg, code);
		stack<reg> stk = strip(interference, n_reg, v_reg);

		// TODO: move this inside strip
		for (reg r = 0; r < v_reg; ++r) {
			if (interference.has(r)) {
				mapping[r].status = color::potential_spill;
				// TODO: maybe remove
				interference.remove(r);
				stk.push(r);
			}
		}

		while (!stk.empty()) {
			const auto s = stk.pop();
			// represents a mask of neighboring registers in use
			reg used = 0;
			for (const auto t : interference.ladj[s]) {
				if (mapping[t].status == color::phys_reg)
					used |= bit(mapping[t].address);
			}
			// equivalent of bitwise NOT when you only consider the n_reg first bits
			const reg free = bits(n_reg) ^ used;
			if (free) {
				mapping[s].status  = color::phys_reg;
				// TODO: use some form of heuristic
				mapping[s].address = some_bit_index(free);
			} else {
				mapping[s].status = color::actual_spill;
				spilled = true;
			}
		}

		std::vector<instr> next_code;
		auto next_v_reg = v_reg;
		enum usage { use, def };
		const auto ref = [&] (reg r, usage u) {
			if (mapping[r].status == color::actual_spill) {
				const reg next_r = next_v_reg++;
				const auto opcode = u == use ? instr::load_local: u == def ? instr::store_local: (__builtin_unreachable(), instr::none);
				next_code.emplace_back(opcode, next_r, r);
				r = next_r;
			}
			return r;
		};
		for (const auto ins : code) {
			reg rd = ins.rd;
			reg rs1 = ins.rs1;
			reg rs2 = ins.rs2;
			instr *patchme;
			// force patchme to be valid during the iteration
			next_code.reserve(next_code.size() + 4);
			switch (ins.opcode) {
			case instr::def:
			case instr::imm:
				patchme = &next_code.emplace_back(ins);
				patchme->rd = ref(rd, def);
				break;
			case instr::req:
				rd = ref(rd, use);
				next_code.emplace_back(ins.opcode, rd);
				break;
			case instr::copy:
			case instr::load:
				rs1 = ref(rs1, use);
				patchme = &next_code.emplace_back(ins.opcode, rd, rs1);
				patchme->rd = ref(rd, def);
				break;
			case instr::add:
				rs1 = ref(rs1, use);
				rs2 = ref(rs2, use);
				patchme = &next_code.emplace_back(ins.opcode, rd, rs1, rs2);
				patchme->rd = ref(rd, def);
				break;
			case instr::store:
				rs1 = ref(rs1, use);
				rd  = ref(rd , use);
				next_code.emplace_back(ins.opcode, rd, rs1);
				break;
			case instr::load_local:
				// maybe could delete the variable if this is executed
				patchme = &next_code.emplace_back(ins.opcode, rd, rs1);
				patchme->rd = ref(rd, def);
				break;
			case instr::store_local:
				rd = ref(rd, use);
				next_code.emplace_back(ins.opcode, rd, rs1);
				break;
			case instr::none:
				assert(false);
			}
		}

		code = next_code;
		v_reg = next_v_reg;
	} while (spilled);
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
	std::cout << v;
	std::cout << m;
	return 0;
}


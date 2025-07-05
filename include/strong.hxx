#pragma once
#include <iostream>

template <
	typename Data,
	typename Tag,
	template <typename> typename ...Mixin>
struct strong_typed : Mixin<strong_typed<Data, Tag, Mixin...>>...
{
	using value_type = Data;
	constexpr explicit strong_typed(value_type const &data) : data(data) {}
	value_type data;
	using base_type = strong_typed;
};

template <typename> struct compare {};
template <typename> struct arith {};
template <typename, typename> struct pointer_arith {};
template <typename> struct print {};

template <typename T>
constexpr auto operator<=>(compare<T> const &l, compare<T> const &r)
{
	return static_cast<T const&>(l).data <=> static_cast<T const&>(r).data;
}

template <typename T>
constexpr auto operator++(arith<T> &i)
{
	auto self = static_cast<T&>(i);
	++self.data;
	return self;
}

template <typename T>
constexpr auto operator+(arith<T> const &l, arith<T> const &r)
{
	return T(static_cast<T const&>(l).data + static_cast<T const&>(r).data);
}

template <typename T>
constexpr auto operator-(arith<T> const &l, arith<T> const &r)
{
	return T(static_cast<T const&>(l).data - static_cast<T const&>(r).data);
}

template <typename T>
constexpr auto operator*(arith<T> const &l, arith<T> const &r)
{
	return T(static_cast<T const&>(l).data * static_cast<T const&>(r).data);
}

template <typename T>
constexpr auto operator/(arith<T> const &l, arith<T> const &r)
{
	return T(static_cast<T const&>(l).data / static_cast<T const&>(r).data);
}

template <typename T>
constexpr auto operator%(arith<T> const &l, arith<T> const &r)
{
	return T(static_cast<T const&>(l).data % static_cast<T const&>(r).data);
}

template <typename P, typename T>
constexpr auto operator+(P const &p, pointer_arith<T, P> const &i)
{
	return p + static_cast<T const&>(i).data;
}

// TODO: find a way to make pointer_arith<P> be a
// template <T> => current_ptr_ari<T, T>

template <typename T>
std::ostream &operator<<(std::ostream &os, print<T> const &p)
{
	return os << static_cast<T const&>(p).data;
}


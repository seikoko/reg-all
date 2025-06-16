#include <concepts>
#include <utility>
#include <climits>
#include <cassert>


namespace bits {

	namespace impl {
		using std::size_t;

		template <typename Expr, typename Output>
		struct full {
			constexpr auto eval(size_t size, Output dest) const {
				for (size_t i = size_t(0); i < size; ++i, ++dest) {
					*dest = static_cast<Expr const&>(*this).partial(i);
				}
			}
		};

		template <typename Expr>
		struct reduce {
			constexpr auto eval(size_t size) const {
				auto accum = Expr::neutral;
				for (size_t i = size_t(0); i < size; ++i) {
					accum = static_cast<Expr const&>(*this).partial(i, std::move(accum));
				}
				return accum;
			}
		};

	}

	using Unit = std::uint64_t;

	struct lazy : impl::full<lazy, Unit*> {
		Unit *data;

		lazy(Unit *data) : data(data) {}
		constexpr auto partial(size_t idx) const { return data[idx]; }
	};

	struct managed {
		Unit *ptr;
		size_t size;

		static size_t unit_size()
		{
			return CHAR_BIT * sizeof *ptr;
		}

		size_t bit2idx(size_t bit) const
		{
			return bit / unit_size();
		}

		explicit managed(size_t bits)
		{
			bits = (bits + unit_size() - 1) / unit_size() * unit_size();
			size = bit2idx(bits);
			ptr  = new Unit[size]{};
		}

		~managed()
		{
			memset(ptr, 0xbe, size * sizeof *ptr);
			delete[] ptr;
		}

		bits::lazy lazy(size_t bit = 0)
		{
			assert(bit2idx(bit) < size);
			return { ptr + bit2idx(bit) };
		}

		managed &set(size_t bit)
		{
			const auto shift = bit % unit_size();
			ptr[bit2idx(bit)] |= Unit(1) << shift;
			return *this;
		}

		managed &clear(size_t bit)
		{
			const auto shift = bit % unit_size();
			ptr[bit2idx(bit)] &= ~(Unit(1) << shift);
			return *this;
		}

		managed &assign(size_t bit, Unit value)
		{
			const auto idx   = bit2idx(bit);
			const auto shift = bit % unit_size();
			ptr[idx] &= ~(Unit(1) << shift);
			ptr[idx] |=   value   << shift ;
			return *this;
		}

		Unit operator[](size_t bit) const
		{
			const auto shift = bit % unit_size();
			return (ptr[bit2idx(bit)] >> shift) & Unit(1);
		}

	};

	template <typename Expr>
	using full = impl::full<Expr, Unit*>;

	template <typename Expr>
	using reduce = impl::reduce<Expr>;

	template <typename Expr>
	concept expr = std::derived_from<Expr, full<Expr>>;

	template <typename Expr>
	concept reduction = std::derived_from<Expr, reduce<Expr>>
		&& requires (Expr) { Expr::init; };

	template <typename L, typename R>
	struct expr_and : full<expr_and<L, R>> {
		L l;
		R r;

		constexpr expr_and(L l, R r) : l(l), r(r) {}
		constexpr auto partial(size_t idx) const { return l.partial(idx) & r.partial(idx); }
	};
	constexpr auto operator&(expr auto l, expr auto r) { return expr_and(l, r); }

	template <typename L, typename R>
	struct expr_or  : full<expr_or <L, R>> {
		L l;
		R r;

		constexpr expr_or (L l, R r) : l(l), r(r) {}
		constexpr auto partial(size_t idx) const { return l.partial(idx) | r.partial(idx); }
	};
	constexpr auto operator|(expr auto l, expr auto r) { return expr_or (l, r); }

	template <typename L, typename R>
	struct expr_xor : full<expr_xor<L, R>> {
		L l;
		R r;

		constexpr expr_xor(L l, R r) : l(l), r(r) {}
		constexpr auto partial(size_t idx) const { return l.partial(idx) ^ r.partial(idx); }
	};
	constexpr auto operator^(expr auto l, expr auto r) { return expr_xor(l, r); }

	template <typename E>
	struct expr_not : full<expr_not<E>> {
		E e;

		constexpr expr_not(E e) : e(e) {}
		constexpr auto partial(size_t idx) const { return ~e.partial(idx); }
	};
	constexpr auto operator~(expr auto e) { return expr_not(e); }

	template <typename E>
	struct reduce_and : reduce<reduce_and<E>> {
		E e;

		static constexpr Unit neutral = Unit(-1);
		reduce_and(E e) : e(e) {}
		constexpr auto partial(size_t idx, Unit accum) const
		{
			return e.partial(idx) & accum;
		}
	};

	template <typename E>
	struct reduce_or  : reduce<reduce_or <E>> {
		E e;

		static constexpr Unit neutral = Unit(0);
		reduce_or (E e) : e(e) {}
		constexpr auto partial(size_t idx, Unit accum) const
		{
			return e.partial(idx) | accum;
		}
	};

	template <typename E>
	struct count : reduce<count<E>> {
		E e;

		static constexpr Unit neutral = Unit(0);
		count(E e) : e(e) {}
		constexpr auto partial(size_t idx, Unit accum) const
		{
			return __builtin_popcountl(e.partial(idx)) + accum;
		}
	};

}


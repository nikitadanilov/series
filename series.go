//
// This is an implementation of power series, as described in
//
//     FUNCTIONAL PEARLS, Power Series, Power Serious, M. Douglas McIlroy
//
//     (http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.333.3156&rep=rep1&type=pdf)
//
//  

package main

import   "fmt"
import . "math/big"

type Series interface {
	Head() *Rat
	Tail() Series
}

func main() {
	z := &Zero{}
	describe("0", z)
	o := &One{}
	describe("1", o)
	x := &X{}
	describe("x", x)
	describe("6 + 6*x", &Scale{rat(6, 1), &Add{&One{}, &X{}}})
	describe("(x + 1)/7", &Scale{rat(1, 7), &Add{&One{}, &X{}}})
	describe("1/17 + x * (1/17 + x * (1/17 + ...", &Const{rat(1, 17)})
	x2 := &Mul{x, x}
	describe("x^2", x2)
	mx := neg(x)
	describe("-x", mx)
	mo := neg(o)
	describe("-1", mo)
	/* 1 - x */
	xm1 := sub(o, x)
	xx := &Mul{&Add{o, x}, xm1}
	describe("(1+x)*(1-x) == 1 - x^2", xx)
	g := &Div{o, xm1}
	describe("1/(1 - x) = 1 + x + x^2 + ...", g)
	gg := &Div{o, &Mul{xm1, xm1}}
	describe("1/((1 - x)^2) = 1 + 2*x + 3*x^2 + ...", gg)
	fmt.Printf("%s .. %s\n",
		eval(g, rat(1, 8)).FloatString(20), rat(8, 7).FloatString(20))
	describe("(1 - y^2)(x^2) = 1 - x^4", &Subst{xx, x2})
	describe("rev(x) = x", rev(x))
	describe("d/dx(x^2) = 2*x", &Deriv{x2, 0})
	describe("Int(x^2)*dx = x^3/3", &Integral{x2, 0})
	mlnx := &Integral{g, 0}
	describe("Int 1/(1 - x)*dx = -ln(1 - x) = x + x^2/2 + x^3/3 + ...", mlnx)
	/*
	 * Build a cyclic structure: expx = 1 + integral(expx)
	 */
	eint := &Integral{nil, 0}
	expx := &Add{o, eint}
	eint.f = expx
	describe("exp", expx)
	exp1 := eval(expx, rati(1))
	fmt.Printf("%s %s\n", exp1, exp1.FloatString(20))
	sinx := &Integral{nil, 0}
	cosx := &Add{o, neg(&Integral{sinx, 0})}
	sinx.f = cosx
	describe("sin", sinx)
	describe("cos", cosx)
	for i := int64(0); i < 10; i++ {
		fmt.Printf("%s\n", eval(&Add{&Mul{sinx, sinx}, &Mul{cosx, cosx}},
			rati(i)))
	}
	/* rev(sinx) = arcsinx */
	arcsinx := rev(sinx)
	describe("arcsin", arcsinx)
	fmt.Printf("%s\n", eval(arcsinx, rati(1)))
	/* https://www.wolframalpha.com/input/?i=integrate+tan%28x%29+dx */
	describe("Int tan(x) dx", &Integral{&Div{sinx, cosx}, 0})
	describe(" .... ", &Deriv{&Integral{&Div{sinx, cosx}, 0}, 0})
	describe(" .... ", &Div{sinx, cosx})
	/* https://www.wolframalpha.com/input/?i=taylor+series+1+-+1%2F%281%2By%29 */
	describe("rev(1/(1-x) - 1) = (1 - 1/(1+y))",
		rev(&Add{&Div{o, sub(o, x)}, mo}))
	/* https://www.wolframalpha.com/input/?i=taylor+series+ln%281%2Bx%29 */
	lnx := rev(sub(expx, o))
	describe("ln(1+x)", lnx)
	/* https://www.wolframalpha.com/input/?i=integrate+exp%28-x%5E2%29 */
	describe("Int(e^(-x^2))dx", &Integral{&Subst{expx, neg(x2)}, 0})
	describe("sin(arcsin)", &Subst{sinx, arcsinx})
	describe("arcsin(sin)", &Subst{arcsinx, sinx})
	/* https://www.wolframalpha.com/input/?i=taylor+series+arcsin%28arcsin%28x%29%29 */
	describe("arcsin(arcsin)", &Subst{arcsinx, arcsinx})
	/* 3 + 1/10*(3 + 1/10 + ...) = 3 + 1/3 */
	fmt.Printf("%s\n", eval(&Const{rat(3, 1)}, rat(1, 10)))
	/* 9 + 1/10*(9 + 1/10 + ...) = 10 */
	fmt.Printf("%s\n", eval(&Const{rat(9, 1)}, rat(1, 10)))
	atan := rev(&Div{sinx, cosx})
	/* \pi = 16 * atan(1/5) - 4 * atan(1/239) */
	pix := sub(&Scale{rati(16), &Subst{atan, &Scale{rat(1, 5), x}}},
		&Scale{rati(4), &Subst{atan, &Scale{rat(1, 239), x}}})
	describe("pi", pix)
	pi := eval(pix, rati(1))
	fmt.Printf("%s %s\n", pi, pi.FloatString(40))
	kprimes := &Prime{2, 0}
	describe("k-primes", kprimes)
	tprimes := &Map{0, func(n int64) *Rat {
		if isprime(n) {
			return rat(1, n)
		} else {
			return rati(0)
		}
	}}
	describe("t-primes", tprimes)
	for i := int64(-10); i <= int64(+10); i++ {
		for _, j := range []int{10, 100, 1000} {
			kp := ev(kprimes, rat(i, 10), j)
			tp := ev(tprimes, rat(i, 10), j)
			fmt.Printf("%2.2d %6.6d %s %s\n",
				i, j, kp.FloatString(40), tp.FloatString(40))
		}
	}
}

func rat(n int64, d int64) *Rat {
	return new(Rat).SetFrac64(n, d)
}

func rati(n int64) *Rat {
	return rat(n, 1)
}

func ev(s Series, x *Rat, n int) *Rat {
	var t Rat
	var v Rat

	t.SetInt64(1)
	for i := 0; i < n; i++ {
		var m Rat

		m.Mul(&t, s.Head())
		v.Add(&v, &m)
		t.Mul(&t, x)
		s = s.Tail()
	}
	return new(Rat).Set(&v)
}

func eval(s Series, x *Rat) *Rat {
	return ev(s, x, 10)
}

func describe(prefix string, s Series) {
	var i int
	fmt.Printf("%s: [", prefix)
	for i = 0; i < 10; i++ {
		fmt.Printf("%s ", s.Head().RatString())
		s = s.Tail()
	}
	fmt.Printf("... ]\n")
}

type Cons struct {
	head *Rat
	tail Series
}

func (self *Cons) Head() *Rat {
	return self.head
}

func (self *Cons) Tail() Series {
	return self.tail
}

type (
	Zero  struct{}
	One   struct{}
	X     struct{}
	Const struct{ v *Rat }
)

func (self *Zero) Head() *Rat {
	return rati(0)
}

func (self *Zero) Tail() Series {
	return self
}

func (self *One) Head() *Rat {
	return rati(1)
}

func (self *One) Tail() Series {
	return new(Zero)
}

func (self *X) Head() *Rat {
	return rati(0)
}

func (self *X) Tail() Series {
	return new(One)
}

func (self *Const) Head() *Rat {
	return self.v
}

func (self *Const) Tail() Series {
	return self
}

type Pair struct {
	l Series
	r Series
}

type Ind struct {
	f Series
	n int64
}

type (
	Add      Pair
	Scale    Cons
	Mul      Pair
	Div      Pair
	Subst    Pair
	Deriv    Ind
	Integral Ind
)

func (self *Add) Head() *Rat {
	return rati(0).Add(self.l.Head(), self.r.Head())
}

func (self *Add) Tail() Series {
	return &Add{self.l.Tail(), self.r.Tail()}
}

func (self *Scale) Head() *Rat {
	return rati(0).Mul(self.head, self.tail.Head())
}

func (self *Scale) Tail() Series {
	return &Scale{self.head, self.tail.Tail()}
}

func (self *Mul) Head() *Rat {
	return rati(0).Mul(self.l.Head(), self.r.Head())
}

func (self *Mul) Tail() Series {
	return &Add{&Scale{self.l.Head(), self.r.Tail()},
		&Mul{self.l.Tail(), self.r}}
}

func (self *Div) Head() *Rat {
	f := self.l.Head()
	g := self.r.Head()
	if f.Cmp(&Rat{}) == 0 && g.Cmp(&Rat{}) == 0 {
		return (&Div{self.l.Tail(), self.r.Tail()}).Head()
	} else {
		return rati(0).Quo(f, g)
	}
}

func (self *Div) Tail() Series {
	q := rati(0).Quo(self.l.Head(), self.r.Head())
	return &Div{&Add{self.l.Tail(), &Scale{q.Neg(q), self.r.Tail()}}, self.r}
}

func (self *Subst) Head() *Rat {
	if self.r.Head().Cmp(&Rat{}) != 0 {
		panic("free term must be 0.")
	} else {
		return self.l.Head()
	}
}

func (self *Subst) Tail() Series {
	return &Mul{self.r.Tail(), &Subst{self.l.Tail(), self.r}}
}

func rev(s Series) Series {
	if s.Head().Cmp(&Rat{}) != 0 {
		panic("free term must be 0.")
	} else {
		comp := &Subst{s.Tail(), nil}
		rs := &Cons{&Rat{}, &Div{&One{}, comp}}
		comp.r = rs
		return rs
	}
}

func neg(s Series) Series {
	return &Scale{rati(-1), s}
}

func sub(a Series, b Series) Series {
	return &Add{a, neg(b)}
}

func (self *Deriv) Head() *Rat {
	return rati(0).Mul(rati(self.n+1), self.f.Tail().Head())
}

func (self *Deriv) Tail() Series {
	return &Deriv{self.f.Tail(), self.n + 1}
}

func (self *Integral) Head() *Rat {
	if self.n == 0 {
		return &Rat{}
	} else {
		return rati(0).Quo(self.f.Head(), rati(self.n))
	}
}

func (self *Integral) Tail() Series {
	ff := self.f
	if self.n > 0 {
		ff = ff.Tail()
	}
	return &Integral{ff, self.n + 1}
}

type Prime struct {
	p int64
	n int64
}

func (self *Prime) Head() *Rat {
	return rat(1, self.p)
}

func (self *Prime) Tail() Series {
	p := self.p
	n := self.n
	for {
		p++
		n++
		if isprime(p) {
			return &Prime{p, n}
		}
	}
}

func isprime(p int64) bool {
	return new(Int).SetInt64(p).ProbablyPrime(1)
}

type Map struct {
	n   int64
	val func(int64) *Rat
}

func (self *Map) Head() *Rat {
	return self.val(self.n)
}

func (self *Map) Tail() Series {
	return &Map{self.n + 1, self.val}
}

import unittest
from Embedding import *


class TestRational(unittest.TestCase):
	def test_add(self):
		x = Rational(11,8)
		y = Rational(53,200)

		s = x + y
		assert(s == Rational(41,25))

	def test_sub(self):
		x = Rational(11,8)
		y = Rational(53,200)

		d = x-y
		assert (d == Rational(111,100))

	def test_inv(self):
		x = Rational(11,8)

		ix = ~x
		assert(ix == Rational(8,11))

	def test_div(self):
		x = Rational(11,8)
		y = Rational(53,200)

		d = x/y
		assert (d == Rational(275,53))


	def test_continued_frac_convergents(self):
		x = Rational(17,39)

		convs = continued_frac_convergents(x)
		econvs = [Rational(0,1),
				  Rational(1,2),
				  Rational(3,7),
				  Rational(7,16),
				  Rational(17,39)]

		assert(convs == econvs)

		x = Rational(275,53)

		convs = continued_frac_convergents(x)
		econvs = [Rational(5,1),
				  Rational(26,5),
				  Rational(83,16),
				  Rational(275,53)]

		assert (convs == econvs)


	def test_continued_frac_nadic_approx(self):
		x = np.pi
		n = 7
		w = 1

		close = continued_frac_nadic_approx(x,n,w)
		exp = n_adic(3,n,0)

		assert(close == exp)

		w = 6

		close = continued_frac_nadic_approx(x,n,w)
		exp = n_adic(22,n,1)

		assert(close == exp)

		w = 10

		close = continued_frac_nadic_approx(x,n,w)
		exp = n_adic(22,n,1)

		assert (close == exp)

		n = 106
		w = 13

		close = continued_frac_nadic_approx(x,n,w)
		exp = n_adic(333,n,1)

		assert (close == exp)

		n = 10
		w = 100

		close = continued_frac_nadic_approx(x,n,w)
		exp = n_adic(312689,n,5)

		assert (close == exp)


if __name__ == '__main__':
	unittest.main()

import unittest
import coding_challenge


class TestGeneSequencer(unittest.TestCase):

	def test1(self):
		self.assertTrue(coding_challenge.main('test1.txt') == 'ATTAGACCTGCCGGAATAC')

	def test2(self):
		pass


if __name__ == '__main__':
	unittest.main()
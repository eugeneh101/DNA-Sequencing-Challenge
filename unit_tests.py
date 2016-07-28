import unittest
import coding_challenge


class TestGeneSequencer(unittest.TestCase):

	def test1(self):
		self.assertEqual(coding_challenge.main(
			'test_cases/sample_example.txt'), 
			'ATTAGACCTGCCGGAATAC')

	def test2(self):
		self.assertEqual(coding_challenge.main(
			'test_cases/sample_example_shuffled.txt'), 
			'ATTAGACCTGCCGGAATAC')

	def test3(self):
		self.assertEqual(coding_challenge.main(
			'test_cases/cannot_detect_left.txt'), 
			'AAAAACTGATAGATTTAAAAAG')

	def test4(self):
		self.assertEqual(coding_challenge.main(
			'test_cases/cannot_detect_right.txt'), 
			'ATTAGGAAAAAGTCGAAAAA')

	def test5(self):
		self.assertEqual(coding_challenge.main(
			'test_cases/cannot_detect_left_or_right.txt'), 
			'AAAAACTGATAGATTTAAAAAGATTAGGTTTTTGTCGTTTTT')

	def test6(self):
		self.assertEqual(coding_challenge.main(
			'test_cases/cannot_detect_left_or_right_messy.txt'), 
			'AAAAACTGATAGATTTAAAAAGATTAGGTTTTTGTCGTTTTT')


if __name__ == '__main__':
	unittest.main()
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <armadillo>
class TestClass
{
	public:
	int value;
	public:
	TestClass(int x):value(x){}
	TestClass  operator + (const TestClass &test)
	{
		TestClass returnClass(0);
		returnClass.value = value + test.value;
		return TestClass(value + test.value);
	}
	TestClass & operator = (TestClass test)
	{
		value = test.value;
		return *this;	
	}
};
int main()
{
	using namespace boost::numeric::ublas;
	matrix<double> m1(3,3);
	matrix<double> m2(3,3);
	for (unsigned i = 0; i < m1.size1(); i++)
		for (unsigned j = 0; j < m1.size2(); j++)
		{
			m1(i, j) = 3 * i + j;
			m2(i, j) = 6 * i + j;
		}
	matrix<double> m3(3,3);
	m3 = m1 - m2;
	std::cout << m3;
	std::cout << m1;
	std::cout << m2;
	std::cout << prod(m1, m2);
	TestClass class1(2);
	TestClass class2(5);
	TestClass class3(0);
	class3 =(class1 + class2);
	std::cout << class3.value;
}

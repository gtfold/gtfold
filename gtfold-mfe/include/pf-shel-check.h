#include <map>
#include "key.h"
class pf_shel_check {
	public:
	int main();
	void clear();

	int numCount;
	pf_shel_check();
	pf_shel_check(int num);
	int count();
	void add(unsigned char type, int i, int j, bool isNumerator);
	int getNumNumerator();
	int getNumDenominator();
	std::map<Key, int> myMap;
};

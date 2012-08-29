#include "pf-shel-check.h"
#include <stdio.h>

/* Documentation of type codes:
0 = u
1 = up
2 = upm
3 = s1
4 = s2
5 = s3
*/

/*int main(){
	pf_shel_check test = pf_shel_check(7);
	test.add(0, 1, 7, false);
	test.add(1, 1, 3, false);
	printf("%d", (int)test.count());
	return 0;
}*/

pf_shel_check::pf_shel_check() 
{
}
pf_shel_check::pf_shel_check(int len) 
{
	add(0, 1, len, true);
}
void pf_shel_check::clear()
{
       myMap.clear(); 
}

int pf_shel_check::count() {
	//printf("There are %d elements in the numerator.\n There are %d elements in the denominator.\n", getNumNumerator(), getNumDenominator());
	//return getNumNumerator() + getNumDenominator();
	return myMap.size();
}

void pf_shel_check::add(unsigned char type, int i, int j, bool isNumerator) {
	Key key = Key(type, i, j);
	std::map<Key, int>::iterator check = myMap.find(key);
	if(check == myMap.end()) {	
		myMap.insert(std::make_pair(key, isNumerator ? 1 : -1));
		if(!isNumerator) printf("Added to denominator first\n");
	}
	else {
		int value = check -> second;
		value = value + (isNumerator ? 1 : -1);
		myMap.erase(key);
		if(value != 0) myMap.insert(std::make_pair(key, value));
	}
}

int getNumNumerator() {
return 0;
}

int getNumDenominator() {
return 0;
}


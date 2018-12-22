#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <queue>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <string>
#include <cstring>
#include <ctime>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <tuple> 
#include <memory>

using namespace std;

template <typename T>
class PrefixSum1D
{
public:
	PrefixSum1D(int size, int init) {
		m_s = vector<T>(size, init);
	}

	T PartSum(int l, int r)
	{
		if (l > r) swap(l, r);
		if (l == 0) return m_s[r];
		return m_s[r] - m_s[l - 1];
	};

	//#TODO あるデータから累積和を計算する関数

	vector<T> m_s;
};


int main() {
	int N;
	cin >> N;

	return 0;
}
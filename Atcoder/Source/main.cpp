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
#include <random>
#include <iterator>
#include <cassert>
using namespace std;

using LL = long long;
using ULL = unsigned long long;
constexpr int INF = 2147483647;
int dx[4] = { 0,1,0,-1 }, dy[4] = { 1,0,-1,0 };

//�ȉ���C++�̗��K�̂��߂֗̕��v���O�����@����܂�֗�����Ȃ����ǁc
//------------------------------------------�֗��^----------------------------------------------
//�񎟌��̓_
//TODO map,set�̂ق����y����
template <typename T>
struct Point2D
{
	T x;
	T y;

	friend ostream& operator<<(ostream& os, const Point2D<T>& point);

	bool operator < (const Point2D& right) const
	{
		return (this->x == right.x) ? (this->y > right.y) : (this->x > right.x);
	}
};

ostream& operator<<(ostream& os, const Point2D<int>& point)
{
	os << "(" << point.x << ',' << point.y << ")";
	return os;
}

//�����ēǂ݂Â炩�����̂ŃG�C���A�X�e���v���[�g��
template<class T>
using StdAlloc = std::allocator<T>;

//�ȒP�ȂQ�����z��
//T��bool���g���Ƃ��������Ȃ�̂ŁAbitset<1>���g������
template<class T>
class Array2D
{
public:
	struct Array2DAccessor
	{
		LL h;
		LL w;

		Array2DAccessor operator+(const Array2DAccessor& right) const {
			return{ h + right.h, w + right.w };
		}

		const static initializer_list<Array2DAccessor> UpDownLeftRight;
	};

public:
	Array2D(size_t height, size_t width, T init)
	{
		data = vector<vector<T>>(height, vector<T>(width, init));
	}

	//���[�v�@�C���f�b�N�X�Ȃ�
	void ForEach(function<void(T&)> func)
	{
		for (auto& row : data) {
			for (T& column : row) {
				func(column);
			}
		}
	}

	//template<class FUNC>
	//void ForEach(FUNC func)
	//{
	//	for (auto& row : data) {
	//		for (T& column : row) {
	//			func(column);
	//		}
	//	}
	//}

	////���[�v�@�C���f�b�N�X����Areturn��������
	////����Ă݂����ǂ߂�����g���ɂ����c
	//void ForEach(	//function<void(ULL, ULL, T&)> func,
	//				function<bool()> IsReturnAndClosureW = []() {return false; },
	//				function<bool()> IsReturnAndClosureH = []() {return false; }
	//			)
	//{
	//	for (ULL h = 0ULL; h < data.size(); h++) {
	//		if (IsReturnAndClosureW(h, w, data[h][w])) {
	//			continue;
	//		}
	//		for (ULL w = 0ULL; w < data[h].size(); w++) {
	//			if (IsReturnAndClosureH(h, w, data[h][w])) {
	//				break;
	//			}
	//			//func(h, w, data[h][w]);
	//		}
	//	}
	//}

	ULL GetHeight() 
	{
		return data.size();
	}

	ULL GetWidth()
	{
		return data.front().size();
	}

	template <class INDEX>
	vector<T>& operator[](INDEX indx)
	{
		return data[indx];
	}

	bool IsValid(const Array2DAccessor& indx) const
	{
		return 0LL <= indx.h && indx.h < (unsigned)data.size() && 0LL <= indx.w && indx.w < (unsigned)data.front().size();
	}

public:
	vector<vector<T>> data;
};

template<class T>
const initializer_list<typename Array2D<T>::Array2DAccessor> Array2D<T>::Array2DAccessor::UpDownLeftRight = { { 0LL, 1LL },{ 0LL, -1LL },{ -1LL, 0LL },{ 1LL, 0LL } };

//------------------------------------------�֗��^----------------------------------------------

//�[���D��T��
//�ċA�Ŏ��������Stack Overflow�ɂȂ�c�@Atcoder�̃X�^�b�N�T�C�Y�͂����H
template<class Container, class Accessor>
tuple<ULL> DepthFirstSearch(
	Container& field,
	const Accessor& start,
	const initializer_list<Accessor>& searchDir,
	const function<void(const Accessor&, Container&)>& Work,
	const function<bool(const Accessor&, Container&)>& IsSearch
)
{
	tuple<ULL> info = make_tuple(0);	//�T���̏��@�Ă��Ɓ[
	queue<Accessor> taskQueue;
	taskQueue.push(start);

	while (!taskQueue.empty()) {
		const auto task = taskQueue.front();
		taskQueue.pop();

		//�e�����֒T��
		for (const auto& df : searchDir) {
			auto indx = task + df;
			if (!IsSearch(indx, field)) {
				continue;
			}

			Work(indx, field);

			//�T���Ɋւ���e������X�V
			get<0>(info)++;
			taskQueue.push(indx);
		}
	}

	return info;
}

//-----------------------------------------�֗��֐�---------------------------------------------
//�f�o�b�O�p����
//random_device rnd;
//mt19937 mt(0);
//uniform_int_distribution<int> rand10_9(0, 10);
//for (unsigned i = 0; i < N; i++) {
//	numbers[i] = (long long)(rand10_9(mt) - 10);
//	//numbers[i] = 0;
//}

//������������쐬
void MakeSubString(const string& str, vector<string>& subStrList, ULL maxLength) 
{
	for (int i = 0; i < str.size(); i++) {
		for (int j = i; j < str.size() && j - i <= maxLength; j++) {
			subStrList.emplace_back(str.substr(i, j - i + 1));
		}
	}
}

//�R���e�i�Ɋi�[���ꂽ�v�f�����ꂼ�ꂢ�����邩
//�������Ԃ͒x��
template< template<class T, class Allocator = StdAlloc<T>> class Container, class T >
void EachElementNum(const Container<T>& container, map<T, LL>& numberCounts)
{
	for (const auto& data : container) {
		numberCounts[data] = (numberCounts.count(data) == 0) ? 1 : ++numberCounts[data];
	}
}

template< template<class T, class Allocator = StdAlloc<T>> class Container, class T >
void EachElementNum(const Container<T>& container, vector<tuple<T, LL>>& numberCounts)
{
	map<T, LL> eachElementNum;
	EachElementNum(container, eachElementNum);
	for (const auto& info : eachElementNum) {
		numberCounts.emplace_back(make_tuple(info.first, info.second));
	}
}

//�ݐϘa
template< template<class T, class Allocator = StdAlloc<T>> class Container, class T >
void PrefixSum(Container<T>& c)
{
	for(auto now = begin(c) + 1, prev = begin(c); now != end(c); now++, prev++){
		*now += + *prev;
	}
}

//�t�����̗ݐϘa
template< template<class T, class Allocator = StdAlloc<T>> class Container, class T >
void ReversePrefixSum(Container<T>& c)
{
	for (auto now = rbegin(c) + 1, prev = rbegin(c); now != rend(c); now++, prev++) {
		*now += +*prev;
	}
}

//���j�[�N�ȗv�f��
template <class T>
unsigned long long UniqueElementNum(const vector<T>& data)
{
	auto copy = data;

	sort(begin(copy), end(copy));
	copy.erase(unique(begin(copy), end(copy)), end(copy));

	return copy.size();

	//set<unsigned long long>(begin(copy), end(copy)).size();	//�P�s�Ł@�}�N���ɂ��Ă��ǂ�
}

//�����Ԃ͈͓̔����ǂ����@����ς�g���Â炢
template <class T>
bool IsInClose(T x, T min, T max) {
	return (min <= x && x <= max);
}

//�l�������ɂȂ��Ă��邩
template<class T>
bool IsInOrder(const T& value1, const T& value2)
{
	return !(value2 < value1);
}

template<class T>
bool IsInOrder(const T& value1, const T& value2, const T& value3)
{
	return !(value2 < value1) && !(value3 < value2);
}

template<class T>
bool IsInOrder(const T& value1, const T& value2, const T& value3, const T& value4)
{
	return !(value2 < value1) && !(value3 < value2) && !(value4 < value3);
}

//�l���͈͓����ǂ���
template<class T>
inline bool IsInRange(const T& min, const T& value, const T& max)
{
	return IsInOrder(min, value, max);
}

//�ŏ����{�� #TODO template��
long long gcd(long long a, long long b)
{
	long long c;

	if (a < b) {
		a += b; b = a - b; a -= b;
	}

	while (b != 0) {
		c = a % b;
		a = b;
		b = c;
	}

	return a;
}

//�f�����ǂ��� #TODO template��
bool IsPrime(ULL num)
{
	if (num < 2) return false;
	else if (num == 2) return true;
	else if (num % 2 == 0) return false; // �����͂��炩���ߏ���

	double sqrtNum = sqrt(num);
	for (ULL i = 3; i <= sqrtNum; i += 2)
	{
		if (num % i == 0)
		{
			// �f���ł͂Ȃ�
			return false;
		}
	}

	// �f���ł���
	return true;
}

//�f���������@��������isqrt(n)�܂Łj
//�����Ƒ������@������炵���c
void DecomposePrime(ULL n, vector<tuple<ULL, ULL>>& result)
{
	vector<ULL> primeList;

	// ���鐔�̏����l
	ULL a = 2;

	// ��n �� a ( n �� a * a ) �̊ԃ��[�v����
	while (n >= a * a) {
		// a �Ŋ���؂ꂽ��Aa �͑f����
		// �����āA�����鐔�� a �Ŋ���
		// a �Ŋ���؂�Ȃ�������A a �� 1 ����������
		if (n % a == 0) {
			primeList.emplace_back(a);
			n /= a;
		}
		else {
			a++;
		}
	}
	// �Ō�Ɏc���� n �͑f����
	primeList.emplace_back(n);


	//�f���������ꂼ�ꂢ�����邩�𒲂ׂ�
	set<ULL> uniquePrimeList(begin(primeList), end(primeList));
	for (const auto& p : uniquePrimeList) {
		ULL primeNum = count(begin(primeList), end(primeList), p);
		result.emplace_back(make_tuple(p, primeNum));
	}
}

//n����r���o���g�ݍ��킹�̐�
template <typename Type>
Type nCr(Type n, Type r) {
	Type ans = 1;
	for (Type i = n; i > n - r; --i) {
		ans = ans*i;
	}
	for (Type i = 1; i < r + 1; ++i) {
		ans = ans / i;
	}
	return ans;
}

//�V�������ɑΉ��������ꍇ�A�ȉ����Q�l�ɂ���
//�R���e�i�֗v�f��ǉ�
//�v���~�e�B�u�^�Ȃ�
template<class Num, class Container>
void InitContainer(Num num, Container& container) 
{
	container.resize(num);
	for (auto& d : container) {
		cin >> d;
	}
}

//�R���e�i�֗v�f��ǉ�
//Container<Point2D<T>>�œ��ꉻ�@�ǂ݂ɂ����c�@�}�N���ɂ���H
template<class Num, template<class T, class Allocator = StdAlloc<T>> class Container, class T2>
void InitContainer(Num num, Container<Point2D<T2>>& firstContainer)
{
	firstContainer.resize(num);
	for (auto& d : firstContainer) {
		cin >> d.x >> d.y;
	}
}

//�R���e�i�֗v�f��ǉ��@���Ăяo������
template<class Num, class First, class... Rest>
void InitContainer(Num num, First& first, Rest&... rest)
{
	InitContainer(num, first);
	InitContainer(num, rest...);
}

//�v�f��N�Ɗe�v�f���R���e�i�Ɋi�[
//�R���e�i�𕡐�����悤�ɂ��Ă݂��@�e���v���[�g�̗��K
template <class Num, class... Container>
void InitNumAndContainer(Num& num, Container&... container)
{
	cin >> num;
	InitContainer(num, container...);
}

//�R���e�i�̓��e��stream�֏o��
//���͂Ă��Ɓ[��cout�։��s���Ȃ���@�ʂ����Ă��̏������͕֗��Ȃ̂��c�H
//#TODO �������₵����I�[�o�[���[�h���₵����
template <template<class T, class Allocator = StdAlloc<T>> class Container, class T>
void ShowContainer(const Container<T>& container)
{
	ostream_iterator<T> outItr(cout, "\n");
	copy(container.begin(), container.end(), outItr);
	cout << endl;
}

//����͂�����Ɓc
template <template<class T, class Allocator = StdAlloc<T>> class Container, class T>
void ShowContainer(const Container<T>& container, const function<T>& func)
{
	for (const auto& e : container) {
		func(e);
	}
	cout << endl;
}

//���K�@�s�P�ʂŕ������ǂݍ��ނ��߂̃C�f�B�I��
void InputStringForLine()
{
	constexpr size_t size{ 42 };
	vector<std::string> input;

	for (string buffer; getline(std::cin, buffer) && input.size() < size; input.push_back(buffer));
}

//-----------------------------------------�֗��֐�---------------------------------------------

int main()
{
	ULL money;
	cin >> money;

	//���I�v��@
	

	//ULL res = money;
	//for (unsigned i = 0; i < money; i++) {
	//	unsigned count = 0;
	//	unsigned t = i;

	//	while (t > 0) {
	//		count += t % 6;
	//		t /= 6;
	//	}

	//	unsigned t2 = money - i;
	//	while (t2 > 0) {
	//		count += t2 % 9;
	//		t2 /= 9;
	//	}

	//	if (res > count) {
	//		res = count;
	//	}
	//}

	//cout << res << endl;

	//���̃A�C�f�B�A�ł͉����Ȃ��c
	////���炩����100000�܂ł̒l�̃e�[�u�������
	//vector<tuple<unsigned, unsigned>> table;
	//unsigned ninePow = 9;
	//while (ninePow < 100000) {
	//	table.push_back(make_tuple(ninePow, 1));
	//	ninePow *= 9;
	//}

	//unsigned sixPow = 6;
	//while (sixPow < 100000) {
	//	table.push_back(make_tuple(sixPow, 1));
	//	sixPow *= 6;
	//}

	//[&]() {
	//	for (unsigned nine = 9; nine < 100000; nine *= 9) {
	//		for (unsigned six = 6; six < 100000; six *= 6) {
	//			table.push_back(make_tuple(six + nine, 2));
	//		}
	//	}
	//}();

	//table.push_back(make_tuple(1, 1));
	//sort(begin(table), end(table), [&](const tuple<unsigned, unsigned>& a, const tuple<unsigned, unsigned>& b) {
	//	return get<0>(a) < get<0>(b);
	//});

	//for (const auto& e : table){
	//	cout << get<0>(e) << "," << get<1>(e) << endl;
	//}

	////�񕪒T���ň����l��������
	//unsigned count = 0;
	//while (money > get<0>(table[0])) {
	//	auto low = partition_point(begin(table), end(table), [&](const tuple<unsigned, unsigned>& e) {
	//		return get<0>(e) < money;
	//	});
	//	unsigned value = get<0>(*(low - (low != begin(table))));
	//	money -= value;
	//	count += get<1>(*(low - (low != begin(table))));

	//	cout << "money " << money << ", value " << value << ", count " << count << endl;
	//}

	//cout << count << endl;

	return 0;
}


//int main()
//{
//	ULL N;
//	cin >> N;
//
//	vector<ULL> primeList;
//
//	//ULL num = 0;
//
//	//while (primeList.size() < 55) {
//	//	if (IsPrime(num)) {
//	//		primeList.emplace_back(num);
//	//	}
//	//	num++;
//
//	//}
//
//	//for (auto prime : primeList) {
//	//	cout << prime << " ";
//	//}
//
//	//
//	ULL num = 2;
//	while (primeList.size() < N) {
//		if (IsPrime(num)) {
//			primeList.emplace_back(num);
//		}
//
//		num++;
//	}
//
//	map<vector<ULL>, ULL> sumTable;
//
//	auto Check = [&]() {
//		for (ULL one = 0; one < primeList.size(); one++) {
//			vector<ULL> elements(5);
//			elements[0] = primeList[one];
//			ULL sum = 0;
//
//			for (ULL two = one + 1; two < primeList.size(); two++) {
//				elements[1] = primeList[two];
//
//				for (ULL three = two + 1; three < primeList.size(); three++) {
//					elements[2] = primeList[three];
//
//					for (ULL four = three + 1; four < primeList.size(); four++) {
//						elements[3] = primeList[four];
//
//						for (ULL five = four + 1; five < primeList.size(); five++) {
//							elements[4] = primeList[five];
//							if (sumTable.count(elements) == 0) {
//								sumTable[elements] = accumulate(begin(elements), end(elements), 0);
//							}
//
//							if (IsPrime(sumTable[elements])) {
//								return false;
//							}
//						}
//					}
//				}
//			}
//		}
//
//		return true;
//	};
//
//	int debugIndx = 0;
//
//	while (!Check()) {
//		while(!IsPrime(num)) {
//			num++;
//		}
//
//		primeList.back() = num;
//
//		debugIndx++;
//		if (debugIndx % 10000 == 0) {
//			cout << debugIndx << ", " << num << endl;
//		}
//
//		num++;
//	}
//
//	for (auto prime : primeList) {
//		cout << prime << " ";
//	}
//
//	return 0;
//}





////Union-Find Tree
////�l�b�g�ɗ����Ă����̂��R�s�y��������
////#TODO ������Ǝg���Â炢�̂ł��������ɒ�������
//template<class T = LL>
//class UnionFind
//{
//public:
//	std::vector<LL> parent;
//	size_t size = 0;
//
//	void init(const size_t size, const T UNIT_LABEL = 0)
//	{
//		this->size = size;
//		this->UNIT = UNIT_LABEL;
//
//		rank.resize(size, 0);
//		diff_to_parent.resize(size, UNIT_LABEL);
//		parent.resize(size);
//
//		iota(begin(parent), end(parent), 0);
//	}
//
//	//����v�f���O���[�v�ɑ��݂��邩����H
//	//���݂���ꍇ�͍��̗v�f��Ԃ��H�f�[�^�\�����X�V�H
//	//#TODO const�ɂ������悤�ȁc
//	LL find(const LL x)
//	{
//		assert_initialized();
//		if (parent[x] == x) {
//			return x;
//		}
//
//		LL root = find(parent[x]);
//		diff_to_parent[x] += diff_to_parent[parent[x]];
//		parent[x] = root;
//
//		return root;
//	}
//
//	//��̗v�f�������O���[�v�ɑ����邩�ǂ����𔻒�
//	bool is_same(const LL x, const LL y)
//	{
//		assert_initialized();
//		return find(x) == find(y);
//	}
//
//	// false if d is invalid
//	bool unite(const LL x, const LL y, const T d = 0)
//	{
//		assert_initialized();
//
//		T diff_from_root_x_to_root_y = (d + diff_to_root(x)) - diff_to_root(y);
//		LL root_x = find(x);
//		LL root_y = find(y);
//
//		//���������ꍇ
//		if (root_x == root_y) {
//			return d == diff(x, y); // invalid input of d
//		}
//
//		//�����قȂ�ꍇ
//		if (rank[root_x] < rank[root_y]) {
//			std::swap(root_x, root_y);
//			diff_from_root_x_to_root_y = this->UNIT - diff_from_root_x_to_root_y;
//		}
//
//		if (rank[root_x] == rank[root_y]) {
//			++rank[root_x];
//		}
//
//		parent[root_y] = root_x;
//		diff_to_parent[root_y] = diff_from_root_x_to_root_y;
//
//		return true;
//	}
//
//	//
//	T diff(LL x, LL y)
//	{
//		assert(is_same(x, y));
//		return diff_to_root(y) - diff_to_root(x);
//	}
//
//private:
//	inline void assert_initialized()
//	{
//		assert(size > 0);
//	}
//
//	T diff_to_root(LL x)
//	{
//		assert_initialized();
//		find(x); // connect directly to root
//		return diff_to_parent[x];
//	}
//
//	std::vector<LL> rank;
//	std::vector<T> diff_to_parent;
//	T UNIT;
//};
//
//
//int main() 
//{
//	//unsigned N, M;
//	//cin >> N >> M;
//	//
//	////#TODO ������info�������ς����
//
//	//struct Info
//	//{
//	//	unsigned L;
//	//	unsigned R;
//	//	int D;
//	//};
//
//	//vector<Info> infos(M);
//	//for (auto& info : infos) {
//	//	cin >> info.L >> info.R >> info.D;
//	//}
//
//	////(��, ��)�Ƃ����֌W��
//	//for (auto& info : infos) {
//	//	if (info.L < info.R) {
//	//		swap(info.L, info.R);
//	//		info.D = -info.D;
//	//	}
//	//}
//
//	////�\�[�g
//	//sort(begin(infos), end(infos), [](const auto& a, const auto& b) {
//	//	if (a.L > b.L) {
//	//		return true;
//	//	}
//	//	else if(a.L == b.L){
//	//		return a.R >= b.R;
//	//	}
//	//	else {
//	//		return false;
//	//	}
//	//});
//
//	//for (const auto& info : infos) {
//	//	cout << info.L << ", " << info.R << ", " << info.D << endl;
//	//}
//
//
//	////���j�[�N�ȗv�f���𐔂���
//	//map<unsigned, unsigned> unique;
//	//for (const auto& info : infos) {
//	//	unique[info.L] = info.D;
//	//	unique[info.R] = info.D;
//	//}
//
//	//cout << "uniqueNum:" << unique.size() << endl;;
//
//	//vector<unsigned> points(unique.size());
//
//	UnionFind<LL> unionFind;
//	LL N, M;
//	cin >> N >> M;
//
//	unionFind.init(N);
//	bool ok = true;
//
//	for(LL i = 0; i < M; i++){
//		LL L, R, D; // L,R: 1-indexed
//		std::cin >> L >> R >> D;
//
//		--L; --R; // 0-indexed
//
//		if (!unionFind.unite(L, R, D)) {
//			ok = false;
//			break;
//		}
//	}
//
//	cout << (ok ? "Yes" : "No") << endl;
//
//	return 0;
//}




//struct Edge {
//	int distance;
//	int weight;
//
//	bool operator<(const Edge& a) const {
//		return weight < a.weight;
//	}
//};





//unsigned N = 10;
////unsigned N = 1 * 3;	//1�b�ȉ�
////unsigned N = 5 * 1000;	//��4�b
////unsigned N = 1 * 100000;	//
////unsigned N = 2 * 100000;	//
////unsigned N = 3;
//vector<long long> numbers(N);



//vector<long long> numbers;
//numbers.emplace_back(0);
//numbers.emplace_back(-1);
//unsigned N = numbers.size();

//long long debugResult = 0;

////�f�o�b�O�p�𒼉�
//for (unsigned i = 0; i < N; i++) {
//	auto sum = numbers[i];
//	if (sum == 0) {
//		debugResult++;
//	}
//	for (unsigned j = i + 1; j < N; j++) {
//		sum += numbers[j];
//		debugResult += (sum == 0) ? 1 : 0;
//	}
//}



//int main()
//{
//	ULL N, C;
//	cin >> N >> C;
//
//	struct Sushi {
//		ULL x;
//		ULL v;
//	};
//
//	vector<Sushi> sushiSet(N);
//	for (auto& sushi : sushiSet) {
//		cin >> sushi.x >> sushi.v;
//	}
//
//	sort(begin(sushiSet), end(sushiSet), [](const auto& a, const auto& b) {
//		//return a.x / a.v < b.x / b.v;
//		return a.v > b.v;
//	});
//
//	for (auto& sushi : sushiSet) {
//		cout << sushi.x << " " << sushi.v << endl;
//	}
//
//
//
//}


//
//struct A {
//public:
//	A(int init) :value(init) {}
//	int value;
//	void print() { std::cout << "I am A, value is " << value << std::endl; }
//};
//
//struct B {
//public:
//	B(int init) :value(init) {}
//	int value;
//	void print() { std::cout << "I am B, value is " << value << std::endl; }
//};
//
//struct C {
//public:
//	C(int init) :value(init) {}
//	int value;
//	void print() { std::cout << "I am C, value is " << value << std::endl; }
//};
//
//struct f {
//	template<typename T>
//	int operator()(T& object, int result, int n) {
//		object.value = (object.value * n) + result;
//		object.print();
//		return object.value;
//	}
//
//	template<typename T>
//	int operator()(T& object, int n) {
//		object.value *= n;
//		object.print();
//		return object.value;
//	}
//
//	template<typename T>
//	void operator()(T& object) {
//		cout << object.value << endl;
//	}
//};


//template<size_t index, size_t end, bool isEnd = index == end>
//struct forwardExecute;
//
//template<size_t index, size_t end>
//struct forwardExecute<index, end, false>
//{
//	template<typename Tuple, typename F, typename... Args>
//	static void Execute(Tuple&& tuple, F&& f, Args&&... args)
//	{
//		f(std::get<index>(tuple), std::forward<Args>(args)...);
//		forwardExecute<index + 1, end>::Execute(
//			std::forward<Tuple>(tuple),
//			std::forward<F>(f),
//			std::forward<Args>(args)...
//		);
//
//	}
//};
//
//template<size_t index, size_t end>
//struct forwardExecute<index, end, true>
//{
//	template<typename Tuple, class F, typename... Args>
//	static void Execute(Tuple&& tuple, F&& f, Args&&... args)
//	{
//		//end
//	}
//};
//
//template< std::size_t begin, std::size_t end, typename Tuple, typename F, typename... Args>
//void forwardExecuteApply(Tuple&& tuple, F&& f, Args&&... args)
//{
//	forwardExecute<begin, end>::Execute(
//		std::forward<Tuple>(tuple),
//		std::forward<F>(f),
//		std::forward<Args>(args)...
//	);
//}

//{
//	//A a(1);
//	//B b(2);
//	//C c(3);
//
//	//auto tuple = std::make_tuple(a, b, c);
//
//	//forwardExecuteApply<0, 3>(tuple, f(), 2);
//	//forwardExecuteApply<0, 3>(tuple, f());
//
//	//auto test = make_tuple(1, 2, false);
//	//forwardExecuteApply<0, 3>(test, t());
//}
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

//以下はC++の練習のための便利プログラム　あんまり便利じゃないけど…
//------------------------------------------便利型----------------------------------------------
//二次元の点
//TODO map,setのほうが楽かも
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

//長くて読みづらかったのでエイリアステンプレートに
template<class T>
using StdAlloc = std::allocator<T>;

//簡単な２次元配列
//Tにboolを使うとおかしくなるので、bitset<1>を使うこと
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

	//ループ　インデックスなし
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

	////ループ　インデックスあり、return条件あり
	////作ってみたけどめっちゃ使いにくい…
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

//------------------------------------------便利型----------------------------------------------

//深さ優先探索
//再帰で実装するとStack Overflowになる…　Atcoderのスタックサイズはいくつ？
template<class Container, class Accessor>
tuple<ULL> DepthFirstSearch(
	Container& field,
	const Accessor& start,
	const initializer_list<Accessor>& searchDir,
	const function<void(const Accessor&, Container&)>& Work,
	const function<bool(const Accessor&, Container&)>& IsSearch
)
{
	tuple<ULL> info = make_tuple(0);	//探索の情報　てきとー
	queue<Accessor> taskQueue;
	taskQueue.push(start);

	while (!taskQueue.empty()) {
		const auto task = taskQueue.front();
		taskQueue.pop();

		//各方向へ探索
		for (const auto& df : searchDir) {
			auto indx = task + df;
			if (!IsSearch(indx, field)) {
				continue;
			}

			Work(indx, field);

			//探索に関する各種情報を更新
			get<0>(info)++;
			taskQueue.push(indx);
		}
	}

	return info;
}

//-----------------------------------------便利関数---------------------------------------------
//デバッグ用乱数
//random_device rnd;
//mt19937 mt(0);
//uniform_int_distribution<int> rand10_9(0, 10);
//for (unsigned i = 0; i < N; i++) {
//	numbers[i] = (long long)(rand10_9(mt) - 10);
//	//numbers[i] = 0;
//}

//部分文字列を作成
void MakeSubString(const string& str, vector<string>& subStrList, ULL maxLength) 
{
	for (int i = 0; i < str.size(); i++) {
		for (int j = i; j < str.size() && j - i <= maxLength; j++) {
			subStrList.emplace_back(str.substr(i, j - i + 1));
		}
	}
}

//コンテナに格納された要素がそれぞれいくつあるか
//処理時間は遅い
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

//累積和
template< template<class T, class Allocator = StdAlloc<T>> class Container, class T >
void PrefixSum(Container<T>& c)
{
	for(auto now = begin(c) + 1, prev = begin(c); now != end(c); now++, prev++){
		*now += + *prev;
	}
}

//逆向きの累積和
template< template<class T, class Allocator = StdAlloc<T>> class Container, class T >
void ReversePrefixSum(Container<T>& c)
{
	for (auto now = rbegin(c) + 1, prev = rbegin(c); now != rend(c); now++, prev++) {
		*now += +*prev;
	}
}

//ユニークな要素数
template <class T>
unsigned long long UniqueElementNum(const vector<T>& data)
{
	auto copy = data;

	sort(begin(copy), end(copy));
	copy.erase(unique(begin(copy), end(copy)), end(copy));

	return copy.size();

	//set<unsigned long long>(begin(copy), end(copy)).size();	//１行版　マクロにしても良い
}

//ある閉区間の範囲内かどうか　やっぱり使いづらい
template <class T>
bool IsInClose(T x, T min, T max) {
	return (min <= x && x <= max);
}

//値が昇順になっているか
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

//値が範囲内かどうか
template<class T>
inline bool IsInRange(const T& min, const T& value, const T& max)
{
	return IsInOrder(min, value, max);
}

//最小公倍数 #TODO template化
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

//素数かどうか #TODO template化
bool IsPrime(ULL num)
{
	if (num < 2) return false;
	else if (num == 2) return true;
	else if (num % 2 == 0) return false; // 偶数はあらかじめ除く

	double sqrtNum = sqrt(num);
	for (ULL i = 3; i <= sqrtNum; i += 2)
	{
		if (num % i == 0)
		{
			// 素数ではない
			return false;
		}
	}

	// 素数である
	return true;
}

//素因数分解　試し割り（sqrt(n)まで）
//もっと早い方法もあるらしい…
void DecomposePrime(ULL n, vector<tuple<ULL, ULL>>& result)
{
	vector<ULL> primeList;

	// 割る数の初期値
	ULL a = 2;

	// √n ≧ a ( n ≧ a * a ) の間ループ処理
	while (n >= a * a) {
		// a で割り切れたら、a は素因数
		// そして、割られる数を a で割る
		// a で割り切れなかったら、 a を 1 増加させる
		if (n % a == 0) {
			primeList.emplace_back(a);
			n /= a;
		}
		else {
			a++;
		}
	}
	// 最後に残った n は素因数
	primeList.emplace_back(n);


	//素因数がそれぞれいくつあるかを調べる
	set<ULL> uniquePrimeList(begin(primeList), end(primeList));
	for (const auto& p : uniquePrimeList) {
		ULL primeNum = count(begin(primeList), end(primeList), p);
		result.emplace_back(make_tuple(p, primeNum));
	}
}

//n個からr個取り出す組み合わせの数
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

//新しい方に対応したい場合、以下を参考にする
//コンテナへ要素を追加
//プリミティブ型など
template<class Num, class Container>
void InitContainer(Num num, Container& container) 
{
	container.resize(num);
	for (auto& d : container) {
		cin >> d;
	}
}

//コンテナへ要素を追加
//Container<Point2D<T>>で特殊化　読みにくい…　マクロにする？
template<class Num, template<class T, class Allocator = StdAlloc<T>> class Container, class T2>
void InitContainer(Num num, Container<Point2D<T2>>& firstContainer)
{
	firstContainer.resize(num);
	for (auto& d : firstContainer) {
		cin >> d.x >> d.y;
	}
}

//コンテナへ要素を追加　を呼び出す処理
template<class Num, class First, class... Rest>
void InitContainer(Num num, First& first, Rest&... rest)
{
	InitContainer(num, first);
	InitContainer(num, rest...);
}

//要素数Nと各要素をコンテナに格納
//コンテナを複数個取れるようにしてみた　テンプレートの練習
template <class Num, class... Container>
void InitNumAndContainer(Num& num, Container&... container)
{
	cin >> num;
	InitContainer(num, container...);
}

//コンテナの内容をstreamへ出力
//今はてきとーにcoutへ改行しながら　果たしてこの書き方は便利なのか…？
//#TODO 引数増やしたりオーバーロード増やしたり
template <template<class T, class Allocator = StdAlloc<T>> class Container, class T>
void ShowContainer(const Container<T>& container)
{
	ostream_iterator<T> outItr(cout, "\n");
	copy(container.begin(), container.end(), outItr);
	cout << endl;
}

//これはちょっと…
template <template<class T, class Allocator = StdAlloc<T>> class Container, class T>
void ShowContainer(const Container<T>& container, const function<T>& func)
{
	for (const auto& e : container) {
		func(e);
	}
	cout << endl;
}

//練習　行単位で文字列を読み込むためのイディオム
void InputStringForLine()
{
	constexpr size_t size{ 42 };
	vector<std::string> input;

	for (string buffer; getline(std::cin, buffer) && input.size() < size; input.push_back(buffer));
}

//-----------------------------------------便利関数---------------------------------------------

int main()
{
	ULL money;
	cin >> money;

	//動的計画法
	

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

	//このアイディアでは解けない…
	////あらかじめ100000までの値のテーブルを作る
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

	////二分探索で引く値を見つける
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
////ネットに落ちていたのをコピペしただけ
////#TODO ちょっと使いづらいのでいい感じに直したい
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
//	//ある要素がグループに存在するか判定？
//	//存在する場合は根の要素を返す？データ構造を更新？
//	//#TODO constにしたいような…
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
//	//二つの要素が同じグループに属するかどうかを判定
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
//		//根が同じ場合
//		if (root_x == root_y) {
//			return d == diff(x, y); // invalid input of d
//		}
//
//		//根が異なる場合
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
//	////#TODO 乱数でinfoをいっぱい作る
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
//	////(大, 小)という関係に
//	//for (auto& info : infos) {
//	//	if (info.L < info.R) {
//	//		swap(info.L, info.R);
//	//		info.D = -info.D;
//	//	}
//	//}
//
//	////ソート
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
//	////ユニークな要素数を数える
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
////unsigned N = 1 * 3;	//1秒以下
////unsigned N = 5 * 1000;	//約4秒
////unsigned N = 1 * 100000;	//
////unsigned N = 2 * 100000;	//
////unsigned N = 3;
//vector<long long> numbers(N);



//vector<long long> numbers;
//numbers.emplace_back(0);
//numbers.emplace_back(-1);
//unsigned N = numbers.size();

//long long debugResult = 0;

////デバッグ用愚直解
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
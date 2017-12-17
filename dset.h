#if !defined(__UNIONFIND_H)
#define __UNIONFIND_H

#include <vector>
#include <atomic>
#include <iostream>

class DisjointSets {
public:
	DisjointSets(uint32_t size) : mData(size), mRank(size) {
		for (uint32_t i = 0; i < size; ++i) {
			mData[i] = (uint32_t)i;
			mRank[i] = 1;
		}
	}

	uint32_t find(uint32_t id) {
		if (id == mData[id])
			return id;
		uint32_t p = find(mData[id]);
		mData[id] = p;
		return p;
	}

	bool same(uint32_t id1, uint32_t id2) {
		return find(id1) == find(id2);
	}

	uint32_t unite(uint32_t id1, uint32_t id2) {
		int p1 = find(id1);
		int p2 = find(id2);
		if (p1 == p2)
			return p1;
		if (mRank[p1] > mRank[p2]) {
			mRank[p1] += mRank[p2];
			mData[p2] = p1;
			return p1;
		}
		mRank[p2] += mRank[p1];
		mData[p1] = p2;
		return p2;
	}
	std::vector<uint32_t> mData, mRank;
};

#endif /* __UNIONFIND_H */

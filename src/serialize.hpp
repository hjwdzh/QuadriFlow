#ifndef SERIALIZE_H_
#define SERIALIZE_H_

#include <Eigen/Core>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include "adjacent-matrix.hpp"
using namespace Eigen;

inline void Read(FILE* fp, Vector2i& m) {
	fread(&m[0], sizeof(int), 1, fp);
	fread(&m[1], sizeof(int), 1, fp);
}

inline void Save(FILE* fp, Vector2i& m) {
	fwrite(&m[0], sizeof(int), 1, fp);
	fwrite(&m[1], sizeof(int), 1, fp);
}

inline void Save(FILE* fp, MatrixXd& m) {
	int r = m.rows(), c = m.cols();
	fwrite(&r, sizeof(int), 1, fp);
	fwrite(&c, sizeof(int), 1, fp);
	std::vector<double> buffer(r * c);
	for (int i = 0; i < r; ++i) {
		for (int j = 0; j < c; ++j) {
			buffer[i * c + j] = m(i, j);
		}
	}
	fwrite(buffer.data(), sizeof(double), r * c, fp);
}

inline void Read(FILE* fp, MatrixXd& m) {
	int r, c;
	fread(&r, sizeof(int), 1, fp);
	fread(&c, sizeof(int), 1, fp);
	std::vector<double> buffer(r * c);
	fread(buffer.data(), sizeof(double), r * c, fp);
	m.resize(r, c);
	for (int i = 0; i < r; ++i) {
		for (int j = 0; j < c; ++j) {
			m(i, j) = buffer[i * c + j];
		}
	}
}

inline void Save(FILE* fp, MatrixXi& m) {
	int r = m.rows(), c = m.cols();
	fwrite(&r, sizeof(int), 1, fp);
	fwrite(&c, sizeof(int), 1, fp);
	std::vector<int> buffer(r * c);
	for (int i = 0; i < r; ++i) {
		for (int j = 0; j < c; ++j) {
			buffer[i * c + j] = m(i, j);
		}
	}
	fwrite(buffer.data(), sizeof(int), r * c, fp);
}

inline void Read(FILE* fp, MatrixXi& m) {
	int r, c;
	fread(&r, sizeof(int), 1, fp);
	fread(&c, sizeof(int), 1, fp);
	std::vector<int> buffer(r * c);
	fread(buffer.data(), sizeof(int), r * c, fp);
	m.resize(r, c);
	for (int i = 0; i < r; ++i) {
		for (int j = 0; j < c; ++j) {
			m(i, j) = buffer[i * c + j];
		}
	}
}


inline void Save(FILE* fp, VectorXi& m) {
	int r = m.rows();
	fwrite(&r, sizeof(int), 1, fp);
	std::vector<int> buffer(r);
	for (int i = 0; i < r; ++i) {
		buffer[i] = m(i);
	}
	fwrite(buffer.data(), sizeof(int), r, fp);
}

inline void Read(FILE* fp, VectorXi& m) {
	int r;
	fread(&r, sizeof(int), 1, fp);
	std::vector<int> buffer(r);
	fread(buffer.data(), sizeof(int), r, fp);
	m.resize(r);
	for (int i = 0; i < r; ++i) {
		m(i) = buffer[i];
	}
}

inline void Save(FILE* fp, VectorXd& m) {
	int r = m.rows();
	fwrite(&r, sizeof(int), 1, fp);
	std::vector<double> buffer(r);
	for (int i = 0; i < r; ++i) {
		buffer[i] = m(i);
	}
	fwrite(buffer.data(), sizeof(double), r, fp);
}

inline void Read(FILE* fp, VectorXd& m) {
	int r;
	fread(&r, sizeof(int), 1, fp);
	std::vector<double> buffer(r);
	fread(buffer.data(), sizeof(double), r, fp);
	m.resize(r);
	for (int i = 0; i < r; ++i) {
		m(i) = buffer[i];
	}
}

inline void Save(FILE* fp, Link& p) {
	fwrite(&p, sizeof(Link), 1, fp);
}

inline void Read(FILE* fp, Link& p) {
	fread(&p, sizeof(Link), 1, fp);
}

inline void Save(FILE* fp, TaggedLink& p) {
	fwrite(&p, sizeof(TaggedLink), 1, fp);
}

inline void Read(FILE* fp, TaggedLink& p) {
	fread(&p, sizeof(TaggedLink), 1, fp);
}

inline void Save(FILE* fp, double& p) {
	fwrite(&p, sizeof(double), 1, fp);
}

inline void Read(FILE* fp, double& p) {
	fread(&p, sizeof(double), 1, fp);
}

inline void Save(FILE* fp, int& p) {
	fwrite(&p, sizeof(int), 1, fp);
}

inline void Read(FILE* fp, int& p) {
	fread(&p, sizeof(int), 1, fp);
}

template <class T, class F>
inline void Save(FILE* fp, std::pair<T, F>& p) {
	fwrite(&p.first, sizeof(T), 1, fp);
	fwrite(&p.second, sizeof(F), 1, fp);
}

template <class T, class F>
inline void Read(FILE* fp, std::pair<T, F>& p) {
	fread(&p.first, sizeof(T), 1, fp);
	fread(&p.second, sizeof(F), 1, fp);
}

template <class T, class F>
inline void Save(FILE* fp, std::map<T, F>& p) {
	int num = p.size();
	fwrite(&num, sizeof(int), 1, fp);
	for (auto& s : p) {
		fwrite(&s, sizeof(s), 1, fp);
	}
}

template <class T, class F>
inline void Read(FILE* fp, std::map<T, F>& p) {
	int num;
	p.clear();
	fread(&num, sizeof(int), 1, fp);
	for (int i = 0; i < num; ++i) {
		std::pair<T, F> m;
		fread(&m, sizeof(m), 1, fp);
		p.insert(m);
	}
}

template<class T>
void Save(FILE* fp, std::vector<T>& p) {
	int num = p.size();
	fwrite(&num, sizeof(int), 1, fp);
	for (auto& q : p) {
		Save(fp, q);
	}
}

template<class T>
void Read(FILE* fp, std::vector<T>& p) {
	int num;
	fread(&num, sizeof(int), 1, fp);
	p.resize(num);
	for (auto& q : p) {
		Read(fp, q);
	}
}

template<class T>
void Save(FILE* fp, std::set<T>& p) {
	std::vector<T> buffer;
	buffer.insert(buffer.end(), p.begin(), p.end());
	Save(fp, buffer);
}

template<class T>
void Read(FILE* fp, std::set<T>& p) {
	std::vector<T> buffer;
	Read(fp, buffer);
	p.clear();
	for (auto& q : buffer)
		p.insert(q);
}

#endif
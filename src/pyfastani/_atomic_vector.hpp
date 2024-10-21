#ifndef _SAFEVEC_HPP
#define _SAFEVEC_HPP

#include <mutex>
#include <vector>


template <class T>
class atomic_vector: public std::vector<T> {
protected:
    std::mutex mutex;
public:
    atomic_vector(): std::vector<T>(), mutex() {}
    void push_back(const T& val) {
        mutex.lock();
        std::vector<T>::push_back(val);
        mutex.unlock();
    }
};

#endif

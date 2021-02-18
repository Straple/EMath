// eclock
// Developed by Mob

#include<chrono>
#include<iostream>

class eclock {
    typedef std::chrono::steady_clock::time_point time_point;

    time_point begin;

    time_point curTime() const {
        return std::chrono::high_resolution_clock::now();
    }

public:
    eclock() {
        reset();
    }
    void reset() {
        begin = curTime();
    }

    // ms
    double count() const {
        std::chrono::duration<double> time = curTime() - begin;
        return time.count() * 1000;
    }
};
std::ostream& operator <<(std::ostream& output, const eclock& clock) {
    return output << clock.count() << "ms";
}

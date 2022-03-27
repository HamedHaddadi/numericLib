# ifndef TIMER_H_
# define TIMER_H_

class Timer {
    public:
        using timerType = std::chrono::high_resolution_clock;
        Timer(const char* funcName):funcName_{funcName}, start_{timerType::now()} {};
        Timer(const Timer &) = delete;
        auto operator=(const Timer &) = delete;
        Timer(Timer&&) = delete;
        auto operator=(Timer &&) = delete;

        ~Timer() {
            auto stop = timerType::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start_).count();
            std::cout <<"ran "<<funcName_<<" in "<<duration<<" microseconds "<<std::endl;
        }

    private:
        const char* funcName_{};
        const timerType::time_point start_{};
};



# endif 
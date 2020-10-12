#pragma once

#include <mutex>
#include <vector>


class OMPExceptionHandler
{
private:
    std::vector<std::exception_ptr> m_exceptions;
    std::mutex m_lock;
public:
    OMPExceptionHandler(){}

    void Rethrow() {
        if (m_exceptions.size() > 0) {
            std::rethrow_exception(m_exceptions.at(0));
        }
    }

    template <typename Function, typename... Parameters>
    void Run(Function f, Parameters... params)
    {
        try {
            f(params...);
        } catch (...) {
            CaptureException();
        }
    }

    void CaptureException() 
    {
        std::unique_lock<std::mutex> guard{ m_lock };
        m_exceptions.push_back(std::current_exception());
    }
};
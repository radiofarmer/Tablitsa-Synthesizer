#pragma once

#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <future>

class VoiceThreadPool
{
public:
  using ThreadTask = std::function<void()>;

  VoiceThreadPool(std::size_t numThreads)
  {
    start(numThreads);
  }

  void PushTask(ThreadTask task)
  {
    {
      std::unique_lock<std::mutex> lock{ mEventMutex };
      mTasks.emplace(std::move(task));
    }

    mEvent.notify_one();
  }

  void JoinAll()
  {
    stop();
  }

  ~VoiceThreadPool()
  {
    stop();
  }

private:
  void start(std::size_t numThreads)
  {
    for (auto i = 0u; i < numThreads; ++i)
    {
      mThreads.emplace_back([=]() {
          while (true) {
            ThreadTask task;

            {
              std::unique_lock<std::mutex> lock{ mEventMutex };

              mEvent.wait(lock, [=]() { return mStop || !mTasks.empty(); });

              if (mStop && mTasks.empty())
                break;

              task = std::move(mTasks.front());
              mTasks.pop();
            }

            task();
          }
        });
    }
  }

  void stop() noexcept
  {
    {
      std::unique_lock<std::mutex> lock{ mEventMutex };
      mStop = true;
    }

    mEvent.notify_all();

    for (auto& thread : mThreads)
      thread.join();
  }

  std::vector<std::thread> mThreads;
  std::condition_variable mEvent;
  std::mutex mEventMutex;

  std::queue<ThreadTask> mTasks;

  bool mStop = false;
};
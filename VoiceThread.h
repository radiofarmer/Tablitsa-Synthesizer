#pragma once

#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>

#include "SynthVoice.h"

class VoiceThreadPool
{
public:
  using ThreadTask = std::function<void(void)>;

  VoiceThreadPool(std::size_t numThreads = std::thread::hardware_concurrency()) : mNumThreads(numThreads)
  {
    Start();
  }

  void AddTask(ThreadTask newTask);

  void Start()
  {
    std::unique_lock<std::mutex> lk(mQueueMutex);
    for (int i{ 0 }; i < mNumThreads; ++i)
    {
      mThreads.push_back(std::thread(std::bind(&VoiceThreadPool::thread_loop, this)));
    }
  }

  void FinishBlock();
  void Stop();

  ~VoiceThreadPool()
  {
    Stop();
  }

private:

  inline void VoiceThreadPool::thread_loop()
  {
    while (!mStop.load())
    {
      ThreadTask task;

      {
        std::unique_lock<std::mutex> lk(mQueueMutex);
        mBlockCV.wait(lk, [this]() { return !mTaskQueue.empty(); });
        task = mTaskQueue.front();
        mTaskQueue.pop();
      }

      task();

      {
        std::unique_lock<std::mutex> lk(mQueueMutex);
        mActiveVoices--;
        mBlockCV.notify_all();
      }
    }
  }

protected:
  const int mNumThreads;
  int mMaxVoices; // TODO: Replace the hardcoded number in the ConditionalVariable wait statement with this variable. Would also be a good idea to make a TablitsaVoiceAllocator derived class.

  // Conditional Variable flags
  std::atomic<bool> mStop{ false };
  bool mBlockFinished{ false };
  int mActiveVoices{ 0 };

  std::mutex mQueueMutex;
  std::condition_variable mBlockCV;

  std::queue<ThreadTask> mTaskQueue;
  std::vector<std::thread> mThreads;

  ThreadTask mExitTask = [](void) {return;};
};
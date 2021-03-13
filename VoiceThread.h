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
  void thread_loop();

protected:
  const int mNumThreads;

  // Conditional Variable flags
  std::atomic<bool> mStop{ false };
  bool mBlockFinished{ false };
  int mActiveVoices;

  std::mutex mQueueMutex;
  std::condition_variable mBlockCV;

  std::queue<ThreadTask> mTaskQueue;
  std::vector<std::thread> mThreads;
};
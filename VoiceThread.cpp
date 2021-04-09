#include "VoiceThread.h"

void VoiceThreadPool::AddTask(ThreadTask newTask)
{
  {
    std::unique_lock<std::mutex> lk(mQueueMutex);
    mTaskQueue.push(newTask);
    mActiveVoices++;
  }
  mBlockCV.notify_one();
}

void VoiceThreadPool::FinishBlock()
{
  std::unique_lock<std::mutex> lk(mQueueMutex);
  mBlockCV.wait(lk, [this]() { return mActiveVoices <= 0 || mStop.load(); });
  mActiveVoices = 0;
  mBlockCV.notify_all();
}

void VoiceThreadPool::thread_loop()
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

void VoiceThreadPool::Stop()
{
  {
    std::unique_lock<std::mutex> lk(mQueueMutex);
    mStop.store(true);

    // Empty queue
    while (!mTaskQueue.empty())
      mTaskQueue.pop();

    for (int i{ 0 }; i < mNumThreads; ++i)
      mTaskQueue.push(mExitTask);

    mBlockCV.notify_all();
  }
  for (std::thread& t : mThreads)
    t.join();
}
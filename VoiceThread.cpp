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
  mBlockCV.notify_all();
  mActiveVoices = 0;
}

void VoiceThreadPool::thread_loop()
{
  while (true)
  {
    ThreadTask task;

    {
      std::unique_lock<std::mutex> lk(mQueueMutex);
      mBlockCV.wait(lk, [this]() { return !mTaskQueue.empty() || mStop.load(); });
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
  std::unique_lock<std::mutex> lk(mQueueMutex);
  mStop.store(true);
  mBlockCV.notify_all();
  for (std::thread& t : mThreads)
    t.join();

  mThreads.clear();
}
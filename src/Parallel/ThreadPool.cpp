#include "ThreadPool.hpp"
#include <chrono>
#include <vector>
#include <utility>   //std::pair


std::vector<std::pair<size_t,size_t>> ThreadPool::split_range(size_t from,
                                                                   size_t to,
                                                                   size_t n_thread)
{   // contains the [from, to) slices
    std::vector<std::pair<size_t,size_t>> coordinate_list(n_thread) ;

    size_t by   = to / n_thread ;
    size_t from_current = from, to_current = 0 ;

    for(size_t n=0; n<n_thread; n++)
    {   from_current = (to_current == 0) ? 0 : (to_current) ;
        to_current   = from_current + by ;
        (to_current + by > to) ? (to_current = to) : 0 ;
        coordinate_list[n] = std::pair<size_t,size_t>(from_current, to_current) ;
    }

    // if there is some remains, distribute it equally over last threads
    size_t remain = to - to_current ;
    if(remain)
    {   size_t from_correction = 0 ;
        size_t to_correction   = 1 ;
        for(size_t n=n_thread-remain; n<n_thread; n++, from_correction++, to_correction++)
        {   coordinate_list[n].first  += from_correction ;
            coordinate_list[n].second += to_correction ;
        }
    }

    return coordinate_list ;
}

ThreadPool::ThreadPool(size_t n_threads, bool debug)
    : queue_task(),
      queue_mutex(),
      queue_open(true),
      debug(debug)
{   assert(n_threads > 0) ;
    this->threads = std::vector<std::thread>(n_threads) ;
    for(size_t i=0; i<n_threads; i++)
    {   this->threads[i] = std::thread(&ThreadPool::thread_routine, this) ; }
}


ThreadPool::~ThreadPool()
{}


size_t ThreadPool::getNThread() const
{   return this->threads.size() ; }


void ThreadPool::addJob(std::function<void()>&& task)
{   // only add a job in the queues if they are open
    if(this->isQueueOpen())
    {   this->lock_mutex_queue() ;
        // this->debug_print(std::string("started adding job")) ;
        this->queue_task.push(task) ;
        // this->debug_print(std::string("ended adding job")) ;
        this->unlock_mutex_queue() ;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(5)) ;
}


void ThreadPool::join()
{   // closes the queues, later call to addJob() will be effect
    this->close_queue() ;

    // joins the threads
    for(auto& thr : this->threads)
    {   if(thr.joinable())
        {   this->debug_print(std::string("joined a thread")) ;
            thr.join() ;
        }
    }
}


void ThreadPool::thread_routine()
{
    this->debug_print(std::string("started")) ;

    while(true)
    {
        // get a function and the arguments value from the queue
        std::function<void()> task ;
        bool has_task      = false ;
        bool is_queue_open = this->isQueueOpen() ;
        std::pair<size_t,size_t> args ;
        this->lock_mutex_queue() ;
        bool is_queue_empty   = this->queue_task.empty() ;
        if(not is_queue_empty)
        {   // this->debug_print(std::string("fetching from queue")) ;
            task = this->queue_task.front() ;
            this->queue_task.pop() ;
            has_task = true ;
        }
        this->unlock_mutex_queue() ;

        // runs the task
        if(has_task)
        {   this->debug_print(std::string("working")) ;
            task() ;
        }
        // exit
        else if(is_queue_empty and not is_queue_open)
        {   break ; }
        std::this_thread::sleep_for(std::chrono::milliseconds(10)) ;
    }
    this->debug_print(std::string("ended")) ;
}


void ThreadPool::lock_mutex_queue()
{   this->queue_mutex.lock() ;
    // this->debug_print(std::string("locked mutex")) ;
}


void ThreadPool::unlock_mutex_queue()
{   this->queue_mutex.unlock() ;
    // this->debug_print(std::string("unlocked mutex")) ;
}


void ThreadPool::open_queue()
{   this->queue_open = true ; }


void ThreadPool::close_queue()
{   this->queue_open = false ; }


bool ThreadPool::isQueueOpen()
{   return this->queue_open ; }


bool ThreadPool::isDebugOn() const
{   return this->debug ; }


void ThreadPool::debug_print(const std::string& msg, std::ostream& out) const
{   if(this->isDebugOn())
    {   std::hash<std::thread::id> hasher ;
        char message[1024] ;
        sprintf(message, "Thread %zu : %s\n",
                hasher(std::this_thread::get_id()), msg.c_str()) ;
        out << message ;
    }
}


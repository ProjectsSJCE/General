    #include <iostream>
    #include <thread>
    using namespace std;
    //This function will be called from a thread

    void call_from_thread(int i, int j, int k) {
        cout<<"I = "<<i<<j<<k;
        std::cout << "Hello, World" << std::endl;
    }

    int main() {
        //Launch a thread
        int i = 10, a= 15, b = 20;
        cout<<"Aaaa";
        std::thread t1(call_from_thread, i, a, b);
        

        //Join the thread with the main thread
        t1.join();
        cout<<"hello\n";
        return 0;
    }

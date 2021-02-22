#include <string>
#include <stack>
#include <chrono>

using namespace std;

struct TimeRecord{
    auto time;
    string name;

    TimeRecord(auto curr_time, string ticName){
        // Init
        time=curr_time;
        name=ticName;
    }
    chrono::duration<num> diff(){
        // Count how long its been since it ticked
        return chrono::duration<num>(chrono::steady_clock::now() - time).count()
    }
}; 

stack<TimeRecord> ticStack;

void tic(string ticName) {
	currentTic.push(new TimeRecord(chrono::steady_clock::now(), ticName));
}
void toc(string ticName){
	cout << "Time to execute " << currentTic.top().name << ": " << currentTic.top().diff() << endl;
	currentTic.pop();
}

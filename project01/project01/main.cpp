#include <iostream>
#include "src/mainapplication.h"

using namespace std;

int main()
{
    MainApplication m;
    m.runApplication();

    return 0;
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope

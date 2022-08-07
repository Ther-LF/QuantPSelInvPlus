#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;
int main(int argc, char const *argv[])
{
    stringstream ss;
    ss<<123<<"+"<<123;
    cout<<ss.str()<<endl;
    return 0;
}

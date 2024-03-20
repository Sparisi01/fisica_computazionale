#include <array>
#include <iostream>
#include <string>

using namespace std;

class Libro {
   private:
    string titolo;

   public:
    Libro(string titolo) {
        this->titolo = titolo;
    }
    string toString() {
        return this->titolo;
    }
};

int main(int argc, char const* argv[]) {
    string test = "ciao";
    Libro l(test);

    cout << test << endl;
    cout << l.toString() << endl;

    test = test.append("wow");

    return 0;
}

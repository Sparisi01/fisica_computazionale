#include <array>
#include <iostream>
#include <string>

using namespace std;

class Libro {
   public:
    int titolo;
    string titolo2;
    Libro(int titolo) {
        this->titolo = 0;
        this->titolo2 = "";
    }
    int toString() {
        return this->titolo;
    }
};

int main(int argc, char const* argv[]) {
    string test = "ciao";
    Libro l(0);

    printf("%d", sizeof(l.titolo));

    printf("Dimensione titolo %d", sizeof(l));

    return 0;
}

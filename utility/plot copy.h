#include <stdio.h>

#include <iostream>
#include <string>

using namespace std;

void executeGNUPlotCommands(string commands[], int nCommands) {
    FILE* gnupipe = NULL;
    gnupipe = _popen("gnuplot -persistent", "w");

    if (gnupipe == NULL) {
        printf_s("Error opening GNUplot pipe");
        return;
    }

    string singleCommandString;
    for (size_t i = 0; i < nCommands; i++) {
        singleCommandString += (commands[i] + "\n");
    }

    cout << singleCommandString;

    // fprintf(gnupipe, "set term qt persist\n");
    // fprintf(gnupipe, "%s\n", gnuCommand);
}
#include <stdio.h>

void executeGNUPlotCommands(const char *gnuPlotCommand) {
    FILE *gnupipe = NULL;
    gnupipe = _popen("gnuplot -persistent", "w");

    if (gnupipe == NULL) {
        printf_s("Error opening GNUplot pipe");
        return;
    }

    // fprintf(gnupipe, "set term qt persist\n");
    fprintf(gnupipe, "%s\n", gnuPlotCommand);
}
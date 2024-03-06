double integraleTrapezi(double a, double b, int N, double f(double)) {
    double h = (b - a) / N;
    double sum = h / 2 * (f(a) + f(b));

    for (size_t i = 1; i < N; i++) {
        sum += h * f(a + i * h);
    }

    return sum;
}

double integraleSimpson(double a, double b, int N, double f(double)) {
    if (N % 2 == 1) N++;
    double h = (b - a) / N;
    double sum = f(a) + f(b);

    for (size_t i = 1; i <= N / 2 - 1; i++) {
        sum += 2 * f(a + 2 * i * h) + 4 * f(a + ((2 * i) - 1) * h);
    }

    return (sum + 4 * f(a + (N - 1) * h)) * h / 3;
}
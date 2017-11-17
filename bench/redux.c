double round256(double x)
{
    double redux = 0x1.8p52 / 256;
    return (x + redux) - redux;
}

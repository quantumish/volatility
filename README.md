# volatility
C++ code that uses the Newton-Raphson method to approximate the implied volatility of a call option. 

## Usage
### Historical
Compile with the following command: 
```sh
g++ volatility.cpp -std=c++17 -o vol
```
Usage (replace `INITIAL_GUESS` with initial volatility guess - recommended value is 1):
```sh
./vol ./cleanestcalls.csv INITIAL_GUESS
```
### Current
To quickly analyze current data, run `get.py` with `python3 get.py`. 

Note that this uses Yahoo Finance's calculated implied volatility 
— previous testing has shown that this program and Yahoo finance produce very similar numbers when given the same input and therefore 
no extra calculation (and data processing) is necessary.

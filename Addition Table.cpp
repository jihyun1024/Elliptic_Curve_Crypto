#include <stdio.h>

// Coefficient of Elliptic Curve
#define A 0b1000

// Point struct
typedef struct _Point {
	int x;
	int y;
}Point;

// Inverse table of each elements in GF(2^4)
// 몇 번째 index에 해당 원소의 역원이 나오는지
int gf_inverse[16] = {
    -1, // 0 has no inverse
    1,  // 1
    9,  // a
    14, // a + 1
    13, // a^2
    11, // a^2 + 1
    7,  // a^2 + a
    6,  // a^2 + a + 1
    15, // a^3
    2,  // a^3 + 1
    12, // a^3 + a
    5,  // a^3 + a + 1
    10, // a^3 + a^2
    4,  // a^3 + a^2 + 1
    3,  // a^3 + a^2 + a
    8   // a^3 + a^2 + a + 1
};

// GF(2^4) multiplication
int gf_mul(int a, int b) {
    int result = 0;
    int modulo = 0b10011;  // Reduction polynomial: x^4 + x + 1

    for (int i = 0; i < 4; i++) {
        if ((b >> i) & 1) {
            result ^= a << i;
        }
    }

    // Modular reduction
    for (int i = 7; i >= 4; i--) {
        if ((result >> i) & 1) {
            result ^= modulo << (i - 4);
        }
    }

    return result & 0b1111;  // keep only 4 bits
}

// Define Point Doubling
Point doubling(Point P1) {
    Point P3;

    if (P1.x == 0) {
        // P_1 * 2 = Point of infinity when P1.x = 0
        P3.x = 0;
        P3.y = 0;
        return P3;
    }
    else {
        int inv = gf_inverse[P1.x];
        int frac = gf_mul(P1.y, inv);
        int lambda = P1.x ^ frac;

        P3.x = (gf_mul(lambda, lambda)) ^ lambda ^ A;
        P3.y = (gf_mul(P1.x, P1.x)) ^ gf_mul(lambda, P3.x) ^ P3.x;
    }

    return P3;
}

// Define Point Addition
Point addition(Point P1, Point P2) {
    Point P3;

    int dx = P1.x ^ P2.x;
    int dy = P1.y ^ P2.y;

    if (P1.x == P2.x) {
        if (P1.y == P2.y) {
            // x, y좌표 모두 같을 때 -> doubling 호출
            return doubling(P1);
        }
        else {
            // x좌표만 같을 때 -> point of infinity
            P3.x = 0;
            P3.y = 0;
            return P3;
        }
    }
    else {
        int dx_inv = gf_inverse[dx];
        int lambda = gf_mul(dy, dx_inv);

        P3.x = (gf_mul(lambda, lambda)) ^ lambda ^ P1.x ^ P2.x ^ A;
        P3.y = (gf_mul(lambda, P1.x ^ P3.x)) ^ P3.x ^ P1.y;
    }

    return P3;
}

// Define Scalar Multiplication (Left-to-Right method)
Point multiplication(Point P, int d) {
    Point Result = P;

    while (d > 0) {
        Result = doubling(Result);
        
        // If d_i == 1, operate addition
        if ((d & 1) == 1) {
            Result = addition(Result, P);
        }
        
        // Right bit shift
        d = d >> 1;
    }

    return Result;
}

// Print as binary
void print_bin(int num) {
	for (int cnt_i = 3; cnt_i >= 0; cnt_i--) {
		printf("%d", (num >> cnt_i) & 1);
	}
}

// Print as point coordinates
void print_point(Point point) {
    printf("(");
    print_bin(point.x);
    printf(", ");
    print_bin(point.y);
    printf(")");
}

// Print table 
void print_table(Point points[], int size)
{
	printf("[Addition Table]\n");
    Point result;

	for (int cnt_i = 0; cnt_i < size; cnt_i++) {
		for (int cnt_j = 0; cnt_j < size; cnt_j++) {
            result = addition(points[cnt_i], points[cnt_j]);

            print_point(points[cnt_i]);
            printf(" + ");
            print_point(points[cnt_j]);
            printf(" = ");
            print_point(result);
            printf("\n");
		}
	}
}

// Check E(F_2^4) is cyclic
void check_cyclic(Point alpha) {
    Point original = alpha;

    printf("[0 multiplication] = point of infinity\n");
    printf("[%d multiplication] = ", 1);
    print_point(alpha);

    printf("\n[%d multiplication] = ", 2);
    alpha = doubling(alpha);
    print_point(alpha);

    // 0P, 1P, 2P는 이미 출력했으니 3P부터 출력
    for (int cnt_i = 3; cnt_i < 21; cnt_i++) {
        printf("\n[%d multiplication] = ", cnt_i);
        alpha = addition(alpha, original);
        print_point(alpha);

        // Point of infinity인 경우
        if (alpha.x == 0b0000 && alpha.y == 0b0000) {
            printf(" ====> point of infinity\n\n");
        }
    }
}

int main() {
    Point points[] = 
    {
        {0b0000, 0b1011}, {0b0001, 0b0000}, {0b0001, 0b0001},
        {0b0010, 0b1101}, {0b0010, 0b1111}, {0b0011, 0b1100},
        {0b0011, 0b1111}, {0b0101, 0b0000}, {0b0101, 0b0101},
        {0b0111, 0b1011}, {0b0111, 0b1100}, {0b1000, 0b0001},
        {0b1000, 0b1001}, {0b1001, 0b0110}, {0b1001, 0b1111},
        {0b1011, 0b0010}, {0b1011, 0b1001}, {0b1100, 0b0000},
        {0b1100, 0b1100}, {0b1111, 0b0100}, {0b1111, 0b1011}
    };

	int size = sizeof(points) / sizeof(points[0]);
    print_table(points, size);
    
    Point P = { 0b1000, 0b0001 };
    check_cyclic(P);

	return 0;
}

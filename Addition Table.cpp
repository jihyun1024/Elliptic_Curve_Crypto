#include <stdio.h>

// Coefficient of Elliptic Curve
#define A 0b1000

// Point struct
typedef struct _Point {
	int x;
	int y;
}Point;

// Inverse table of each elements in GF(2^4)
// 몇 번째에 해당 원소의 역원이 나오는지
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
    int modulo = 0b10011;  // x^4 + x + 1

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

// Define Point Addition
Point addition(Point P1, Point P2) {
    Point R;

    if (P1.x == P2.x && P1.y == P2.y) {
        // 두 좌표가 같을 경우는 (0, 0)으로 정의
        R.x = 0;
        R.y = 0;
    }
    else {
        int dx = P1.x ^ P2.x;
        int dy = P1.y ^ P2.y;

        if (dx == 0) {
            R.x = 0;
            R.y = 0;
            return R;
        }
        else {
            int dx_inv = gf_inverse[dx];
            int lambda = gf_mul(dy, dx_inv);

            R.x = gf_mul(lambda, lambda) ^ lambda ^ P1.x ^ P2.x ^ A;
            R.y = gf_mul(lambda, P1.x ^ R.x) ^ R.x ^ P1.y;
        }
    }

    return R;
}

// Print as binary
void print_bin(int num) {
	for (int cnt_i = 3; cnt_i >= 0; cnt_i--) {
		printf("%d", (num >> cnt_i) & 1);
	}
}

// Print table 
void print_table(Point points[], int size)
{
	printf("[Addition Table]\n");

	for (int cnt_i = 0; cnt_i < size; cnt_i++) {
		for (int cnt_j = 0; cnt_j < size; cnt_j++) {
			Point result = addition(points[cnt_i], points[cnt_j]);

			printf("(");
			print_bin(points[cnt_i].x);
			printf(", ");
			print_bin(points[cnt_i].y);
			printf(") + (");
			print_bin(points[cnt_j].x);
			printf(", ");
			print_bin(points[cnt_j].y);
			printf(") = (");
			print_bin(result.x);
			printf(", ");
			print_bin(result.y);
			printf(")\n");
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
	return 0;
}

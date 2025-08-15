//g++ -O3 -std=c++11 -fopenmp VSSort.cpp -o VSSort


#include <tuple>
#include <vector>
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <climits>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <random>
#include <queue>
#include <deque>
#include <list>
#include <map>
#include <array>
#include <chrono>
#include <fstream>
#include <functional>
#include <unordered_map>
#include <omp.h>

using namespace std;

#define N 10000000
#define LIM (1<<30)

int y[N];
int tmp[2][N];
int arr[5][N];
int br[2][N/2];
int y2[2][N/2];

int compare(const void* a, const void* b) {
    int x = *(const int*)a;
    int y = *(const int*)b;
    return (x > y) - (x < y);  // オーバーフローしない
}

int rnd(int mini, int maxi) {
	static mt19937 mt((int)time(0));
	uniform_int_distribution<int> dice(mini, maxi);
	return dice(mt);
}

void GenRandom(int a[],int n,int lim){

    for(int i=0;i<n;i++){
        a[i]=rnd(-lim,lim);
    }
    
}

bool check_sort(int a[],int b[],int n){
    int cnt=0;
    for(int i=0;i<n;i++){
        if(a[i]!=b[i]){cnt++;}  
    }

    if(cnt==0){return true;}
    else{return false;}
}

void r_sort(int a[], int n, int y3[]){
    int maxval = 0;
    for (int i = 0; i < n; i++) {
        if (abs(a[i]) > maxval) maxval = abs(a[i]);
    }
    int maxbit = 0;
    while (maxval > 0) {
        maxbit++;
        maxval >>= 8;
    }
    if (maxbit == 0) maxbit = 1;

    for (int loop = 0; loop < maxbit; loop++) {
        int bucket[512] = {0};

        for (int i = 0; i < n; i++) {
            uint32_t key = ((uint32_t)a[i]) ^ 0x80000000u;
            int index = (key >> (loop * 8)) & 0xFF;
            index += 256;  // オフセット
            bucket[index]++;
        }
        for (int i = 1; i < 512; i++) {
            bucket[i] += bucket[i-1];
        }
        for (int i = n - 1; i >= 0; i--) {
            uint32_t key = ((uint32_t)a[i]) ^ 0x80000000u;
            int index = (key >> (loop * 8)) & 0xFF;
            index += 256;
            y3[--bucket[index]] = a[i];
        }
        memcpy(a, y3, n * sizeof(int));
    }
}

const int INSERTION_SORT_THRESHOLD = 16;

void insertionSort(std::vector<int>& arr, int left, int right) {
    for (int i = left + 1; i <= right; i++) {
        int key = arr[i];
        int j = i - 1;
        while (j >= left && arr[j] > key) {
            arr[j+1] = arr[j];
            j--;
        }
        arr[j+1] = key;
    }
}

// Median-of-three to choose a good pivot.
int medianOfThree(std::vector<int>& arr, int left, int right) {
    int mid = left + ((right - left) >> 1);
    if (arr[right] < arr[left])
        std::swap(arr[left], arr[right]);
    if (arr[mid] < arr[left])
        std::swap(arr[mid], arr[left]);
    if (arr[right] < arr[mid])
        std::swap(arr[right], arr[mid]);
    return mid;
}

// Tuned quicksort implementation.
void quickSort(std::vector<int>& arr, int left, int right) {
    if (right - left <= INSERTION_SORT_THRESHOLD) {
        insertionSort(arr, left, right);
        return;
    }
    
    int mid = medianOfThree(arr, left, right);
    // Place pivot at right-1 position.
    std::swap(arr[mid], arr[right - 1]);
    int pivot = arr[right - 1];
    
    int i = left, j = right - 1;
    while (true) {
        while (arr[++i] < pivot) { }
        while (arr[--j] > pivot) { }
        if (i < j)
            std::swap(arr[i], arr[j]);
        else
            break;
    }
    std::swap(arr[i], arr[right - 1]);
    
    quickSort(arr, left, i - 1);
    quickSort(arr, i + 1, right);
}

void fastIntSort(std::vector<int>& arr) {
    if (!arr.empty())
        quickSort(arr, 0, arr.size() - 1);
}

bool sub(vector<int>arr){
    
    // Sort the array using our fast integer sort.
    fastIntSort(arr);
    
    // Verify that the array is sorted.
    bool sorted = true;
    for (int i = 1; i < N; ++i) {
        if (arr[i-1] > arr[i]) {
            sorted = false;
            break;
        }
    }
    return sorted;
}

int ms[N];

int main() {
    
    double start;
 
    GenRandom(arr[0],N,LIM);
    for(int i=0;i<N;i++){
        arr[1][i]=arr[0][i];
        arr[2][i]=arr[0][i];
        arr[3][i]=arr[0][i];
        arr[4][i]=arr[0][i];
        if(i<N/2){
            br[0][i]=arr[0][i];
        }
        else{
            br[1][i-(N/2)]=arr[0][i];
        }
    }

    vector<int>vec;

    for(int i=0;i<N;i++){
        vec.push_back(arr[0][i]);
    }

    start = omp_get_wtime();

    qsort(arr[0], N, sizeof(int), compare);

    start = omp_get_wtime()-start;

    printf("qsort=%lf\n",start);

    //start = omp_get_wtime();

    //radix_sort(vec.begin(),vec.end());

    //qsort(arr, N, sizeof(int), compare);

    //start = omp_get_wtime()-start;

    //printf("radixsort=%lf\n",start);

    start = omp_get_wtime();

    sort(arr[1],arr[1]+N);

    //qsort(arr, N, sizeof(int), compare);

    start = omp_get_wtime()-start;

    printf("sort=%lf\n",start);

    start = omp_get_wtime();
   
    r_sort(arr[2],N,y);

    start = omp_get_wtime()-start;

    printf("radix_sort=%lf\n",start);

    start = omp_get_wtime();

    bool ok=sub(vec);
    //r_sort(narr,N);

    start = omp_get_wtime()-start;

    printf("quick_sort=%lf\n",start);

    start=omp_get_wtime();

    #pragma omp parallel for
    for(int i=0;i<2;i++){
    r_sort(br[i],N/2,y2[i]);
    }
    
    merge(br[0], br[0] + (N/2), br[1], br[1] + (N/2), ms);

    start = omp_get_wtime()-start;
    
    printf("MT_radix_sort=%lf\n",start);

    bool j[3]={false,false,false};

    start= omp_get_wtime();

    j[0]=check_sort(arr[0],arr[1],N);    

    start = omp_get_wtime()-start;

    printf("O(N)=%lf\n",start);

    j[1]=check_sort(arr[0],arr[2],N);
    j[2]=check_sort(arr[1],arr[2],N);
    bool MT=check_sort(ms,arr[0],N);

    if(j[0]&&j[1]&&j[2]&&ok&&MT){printf("ok\n");}
    else{printf("ng\n");}
    
    return 0;
}

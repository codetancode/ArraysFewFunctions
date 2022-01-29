import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Set;
import java.util.List;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;


interface ArraySol<T>{
    public void showArray(T[] arr); // T_O(n),
    public void show2dArr(int[][] arr);
    public void show(int[] arr);//T_O(n), 
    public T[] reverse(int start, T[] arr, int end);//T_O(n/2)=>O(n)
    public int[] maxminPair(Integer[] arr);
    public T kthmax(T[] arr, int k);//T_O(nlogn)//arraydefault sort
    public T kthmin(T[] arr, int k);//T_O(nlogn)//arraydefault sort
    //sort 3 kind of values in an array
    public int[] sort1s2s3s(int[] arr);//T_O(2n)=>O(n) 
    public void neg1side(int[] arr);//T_O(n)//one extra variable as pointer
    public int unionCountOfArray(T[] a, T[] b);//using set T_O(n)
    public T[] unionOfArray(T[] a, T[] b);//using set T_O(n)
    public int intersectCountOfArray(int[] a, int[] b);//using set T_O(n)
   //public int[] intersectOfArray(int[] a, int[] b);//using set T_O(n)
    public void rotate1Arrayleft(int[] arr);//T_O(n)
    public void rotate1Arrayright(int[] arr);//T_O(n)
    //kadanesAlgo
    public int maxSubarraySum(int[] arr);//using standard algo T_O(n)
    //window of size(sliding window)
    public int maxKSubarraySum(int[] arr, int size);//using sliding window T_O(n)
    public int minKSubarraySum(int[] arr, int size);//using sliding window T_O(n)
    //minimize height of array(min height difference(max-min) after +/- k in every element)
    public int minimizedHeight(int[] arr, int k);//using sort, then max min default function T_O(nlogn + n-2)
    //minimum jumps to reach the end of array
    //traverse and keep max variable scope , reduce steps id all steps are done take jump++
    public int minJumpsToEnd(int[] arr);//in T_O(n),
    //find duplicate in distinct value array (only 1 value is twice)
    public int duplicate2Val(int[] arr);//T_O(n) use set
    //merge and do inplace chages merge without extra memory(so that arr1 has smallest number and arr2 has bigger numbers in assending order)
    public void mergeInplace(int[] arr1, int[] arr2);//T_O(nlogn + mlogm)
    //merge intervals in 2d array
    //esential code to sort2d array based on 1st element in internal array
    //Arrays.sort(intervals, (a, b) -> a[0]-b[0]); #==> comparator function of (a, b) return difference of a[0]-b[0] where a and b are intervals arrays
    //similarly for decreasing sort 
    //Arrays.sort(intervals, (a, b) -> b[0] - a[0])
    //using getLast() function of LinkedList helps
    //converting LinkedList to array of primary type 
    //list.toArray(new int[list.size()][]);# to convert list of int[] into int[][]
    public int[][] mergeIntervals(int[][] intervals);//T_O(n) 
    //next permutation in lexological order if heighest permutation i.e decending order sort in ascending order
    //goal is to find 1st decresing element from right , 
    //then find bigger number than that from right, 
    //swap them and reverse the array from new biggers newxt index to end
    //else if no smallest element found from right i.e array in decending order, so just return reverse(sorted ascending)
    public void nextPermutation(int[] arr);
    //count inversion .. i.e no. of swap for a sorted array(consider all elements are distinct)
    //use twaked verion od merge sort and return count if left element is greated that right
    public int inversionCount(int[] arr);//T_O(nlogn) space O(n)
    //can do in 2 ways using hash map and taking a copy of sorted array and traversing array maping indices with sorted array and updating the hashmap
    //can use tweak merge sort , to return count if left element is greater than right element 
    //best time to buy stocks consider all positive value(no going back)
    public int bestProfit(int[] prices);//T_O(n)
    //traverse keeping track of min_cost and max profit variable
    //count inversion-swaps required to convert array in sorted array, 
    //get pair count
    public int countPairSum(int[] arr, int sum);//T_O(n+n)=>T_O(n)
    //can use hashmap to get the frequency first then check by traversing if(sum-element) is there in hashmap or not
    //if there count += hashmap.get(sum - arr[i]) , add count the frequency of it in array(as we are not taking ount distinct pairs)
    //this will the count of all the elements whoes sum could be sum , so return count/2, if number gets repeated i.e(sum-arr[i] == arr[i])
    //reduce the count by one 
    //get maximun product sub array
    //subArrayHasSum0 return true or false
    public boolean subArrayHasSum0(int[] arr);//O(n) 
    public List<Integer> commonElementIn3Sorted(int[] a,int[] b, int[] c);//T_O(n+m+l)
    //3 pointers incrementing if pointing to less value
    //rearange array in alternating +/- values
    public void rearrangeAlternate(int[] arr);
    //take 3 pointers 1 at odd position for -ve swap, 2 at even position making sure positive or swap
    // 3rd to iterate through the array ane swaping//mind the length while doing +=2 in 1st 2 pointer
    //consider 500 size int array and do the calculation , ex output for fact(100)//try doing it with char array to reduce memory
    public void longFactorial(int n);
    //take a result array of 500 size, keep multipling from 2 to =n and storing in result
    //mind carry in loop while multipling each time from 0 to size of result  
    //longest consecutive subsequence(array is un sorted)--try to do in O(n) using hashSet, to find pattern from (!set.contains(element-1)) 
    //to where ever it goes actual length is (incrementation) - start point
    public int lenLongConsecutiveSub(int[] arr);
    //sort array in nlogn then count consequtive max(remember could have duplicate in array)
    //max sub array product
    public int maxSubProd(int[] arr);
    //could be -ve elements
    //use 3 var ,max min and actual max while traversing(math.max, Math.min)
    //frquency of elements more that n/k in an array
    public void freqKMore(int[] arr, int k);
    //get all elements having frequency more that n/k times in an array
    //there are many approches, sort then count in O(nlogn)
    //use hashmap in O(n)
    // try to do it in O(nk). using tetrics way having array of k-1 elements(storing elements and count)
    //iterate throught the array add/remove count++/count--(all), remove when no place the for new element, before adding 
    //new element reduce count of all by 1, after done with the array count numbers in array
    public int maxProfitAtAll(int[] prices);
    //add positive as profit (2nd number-1st) in traversal
    //use peek vally way       8  
    //                    5   /
    //                  3   4
    //#public boolean isASubsetOf(int[] a, int[] b);
    //find the min length check if all of small numb in other array
    public boolean ifTripletsum(int[] arr, int sum);//T_(n2)/M_(1)without hash set, can use sort
    //in 2 loops use 3 pointer 2 at start and 1 at last
    //trapping rain water return ware traped in the heights[2, 0, 3, 0, 4]
    public int trapingRain(int[] arr);//T_O(n)/M_O(n)
    // example array [2, 0, 1, 0, 3]-> water trapped is 5 unit
    //we ned to know max left height, max right then water is min of (left , right)height - current height
    // main formula for that is water+= min(leftheight, right height) - current height
    //idea is to have left array and right array stor the max height in the indices, then traver through array for formula
    //public int countMinSubForSumGreaterThanX(int[] arr, int x);//T_O(n)
    //idea in a loop till last element , have 2 pointer increase the end till sum>x, then while sum>x reduce st
    //public void partationInto3(int[] arr, int a, int b);//T_O(n) where a<b
    //idea traverse from st to last if <a swap increment i, if > b swap from last decrese last, else continue iterating in the array
    //public minSwapKtogether(int[] arr, int k);T_O(n)
    //idea count all <=k elements, make window of that, count >k in the window, keep min of the that in final variable, traverse the array with(count window)
    //changing reqSwaps . return min of req swaps 
    //public void spiralTraversal(int[][] matrix);//for nxm size 
    public int[] solve(int[] A, int[] B);

}
class ArraySolClass<T> implements ArraySol<T>{
        
    public void showArray(T[] arr){
        System.out.print("\n[ ");
        for (int i = 0;i < arr.length - 1; i++ ){
            System.out.print(arr[i]+", ");
        }
        System.out.print(arr[arr.length - 1]+ " ]" + "(" + arr.length + ")");
    }

    public void show(int[] arr){
        System.out.print("\n[ ");
        for (int i = 0;i < arr.length - 1; i++ ){
            System.out.print(arr[i]+", ");
        }
        System.out.print(arr[arr.length - 1]+ " ]" + "(" + arr.length + ")");
    }

    public void show2dArr(int[][] arr){
        System.out.print("\n[ ");
        for(int i = 0;i < arr.length;i++){
            for (int j = 0;j < arr[i].length;j++){
                System.out.println(arr[i][j]);
            }
            System.out.println(" ]");
        }
        
    } 
   
    public T[] reverse(int start, T[] arr, int end){
        T temp;
        while( start < end ){
            temp = arr[start];
            arr[start] = arr[end];
            arr[end] = temp;
            start++;
            end--;
        }
        return arr;
    } 

    public void swap(int[] arr, int i, int j){
        int temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }

    public int[] maxminPair(Integer[] arr){
        int[] maxmin = new int[2]; 
        if ( arr.length == 1){
            maxmin[0] = arr[0];
            maxmin[1] = arr[0];
            return maxmin;
        }
        //initiate the array with some numbers
        if(arr[0] < arr[1]){
            maxmin[0] = arr[1];
            maxmin[1] = arr[0];
        }
        else{
            maxmin[0] = arr[0];
            maxmin[1] = arr[1];
        }
        //start from 2 as already initiated
        for (int i = 2;i < arr.length; i++){
            if (maxmin[0] < arr[i]){
                maxmin[0] = arr[i];
            }
            else if (maxmin[1] > arr[i]){
                maxmin[1] = arr[i];
            }
        }
            return maxmin;    
        }
        
    public T kthmax(T[] arr, int k){
        //sort then pick index
        //arrays.sort in java used dualt pivotal quick sort so expected(nlogn)
        Arrays.sort(arr, Collections.reverseOrder());
        return arr[k-1];
 
    }

    public T kthmin(T[] arr, int k){
        //sort then pick index
        //arrays.sort in java used dualt pivotal quick sort so expected(nlogn)
        Arrays.sort(arr);
        return arr[k-1];
 
    }

    public int[] sort1s2s3s(int[] arr){
        int c0 = 0;
        int c1 = 0;
        int c2 = 0;
        int i = 0;
        for(i = 0;i < arr.length; i++){
         if (arr[i] == 0){
             c0++;
         }
         else if (arr[i] == 1){
             c1++;
         }   
         else if (arr[i] == 2){
             c2++;
         }
         else{
             System.out.println("Array not only have 0s, 1s, 2s");
         }
        }
            i = 0;
            while(c0 > 0){
                arr[i++] = 0;
                c0--;
            }
            while(c1 > 0){
                arr[i++] = 1;
                c1--;
            }
            while(c2 > 0){
                arr[i++] = 2;
                c2--;
            }
        
        return arr;
        }

    public void neg1side(int[] arr){
        if (arr.length > 0){
            int neg = 0;
            int temp;
            //swaping ig -ve encountered and keeping track if neg index with neg
            for (int i = 0;i < arr.length; i++){
                if (arr[i] < 0){
                    temp = arr[i];
                    arr[i] = arr[neg]; 
                    arr[neg] = temp;
                    neg++;
                    
                }
            }

        }
    }



    
    public int intersectCountOfArray(int[] a, int[] b){
        //array can have duplicates and be unsorted
        //sorting array and finding common elements count
        Arrays.sort(a);
        Arrays.sort(b);
        int com = 0;
        int i = 0;
        int j = 0;
        while(i< a.length && j < b.length){
            if (a[i] > b[j]){
                j++;
            }
            else if(a[i] < b[j]){
                i++;
            }
            else{
                //if a[i]==a[j]
                //checking if issue of -ve index refference to check sequence
                if (i ==0 && j==0){
                    com++;
                }
                else{
                    if (i > j){
                        //can check i-1 for sequence
                        if (a[i] != a[i-1]){com++;}
                    }
                    else if(i < j){ 
                        //can check i-1 for sequence
                        if (b[j] != b[j-1] ){com++;}
                }
                else{
                    //i = j but > 0 so can check with any of them for sequence
                    if ((a[i] != a[i-1]) || (b[j] != b[j-1])){com++;}
                }
                }
                
                i++;
                j++;
            }
        }
        return com;

    }

    public int[] solve(int[] A, int[] B) {
        int sum = A[0];
        int i = 0;
        int j = 0;
        int[] c = new int[B.length];
        boolean min = false;
        while(i < B.length){
            while(sum < B[i]){
            	if(j+1 < A.length) {
            		sum += A[++j];
            		System.out.print("Increasing..i sum j "+i+" "+sum+" "+j+"\n");
            	}
            	else {
            		c[i] = A.length;
            		System.out.print("Max out breaking after i sum j "+i+" "+sum+" "+j+"\n");
            		break;
            	}
            	
            }
             
            while(sum > B[i]){
            	if(j-1 >= 0) {
            		sum -= A[j--];
            		System.out.print("Decreasing..i sum j "+i+" "+sum+" "+j+"\n");
            	}
            	else {
            		c[i] = 0;
            		min = true;
            		System.out.print("Min out breaking after i sum j "+i+" "+sum+" "+j+"\n");
            		break;
            	}
            	
                
            }
            this.show(c);
           
            if(c[i] == 4 || !min) {
            	c[i] = j+1;
            	
            	this.show(c);
            }
            i++;
            min = false;
        }
        return (c);
    }
    
    
    public int unionCountOfArray(T[] a, T[] b){
        Set<T> set = new HashSet<T>();
        for (int i = 0;i < a.length ; i++){
            set.add(a[i]);
        }
        for (int i = 0;i < b.length ; i++){
            set.add(b[i]);
        }
        return set.size();

    }

    public T[] unionOfArray(T[] a, T[] b){
        Set<T> set = new HashSet<T>();
        for (int i = 0;i < a.length ; i++){
            set.add(a[i]);
        }
        for (int i = 0;i < b.length ; i++){
            set.add(b[i]);
        }
        return set.toArray(a);
        
    }

    public void rotate1Arrayleft(int[] arr){
        int temp = arr[0], i;
        for (i = 0;i < arr.length - 1;i++){
            arr[i] = arr[i + 1];
        }
        arr[arr.length - 1] = temp;
    }

    public void rotate1Arrayright(int[] arr){
        int temp = arr[arr.length - 1], i;
        for (i = arr.length - 1;i > 0;i--){
            arr[i] = arr[i - 1];
        }
        arr[0] = temp;
    }

    //kadanes algo
    public int maxSubarraySum(int[] arr){
        int sum = 0;
        int max = 0;
        for (int i = 0;i < arr.length; i++){
            sum += arr[i];
            if (max < sum){
                max = sum;
            }
            //if due to -ve integers in between array, put sum = 0,
            //su that new sub sequence can be checked, as older is already in max
            if( sum < 0){
                sum = 0;
            }
        }
        return max;
    }

    public int maxKSubarraySum(int[] arr, int size){
        int max = 0;
        if (arr.length < size || arr.length == 0){
            return -1;
        }
        //if size is same as array len then sum of array is max/min
        if (arr.length == size){
            for (int i = 0;i < size;i++){
                max += arr[i];
            }
        }
        else{
            int st = 0;
            int end = 1;
            int sum = arr[st];
            //get window size sum first then will traverse array 
            while(end < size){
                sum += arr[end++];
            }
            max = sum;
            //now traverse with window on  array 
            while(end < arr.length){
                //changing sum value with incrementing st and end at the same time
                sum -= arr[st++];
                sum += arr[end++];
                if( max < sum){
                    max = sum;
                }
            }

        }

        return max;
    }

    public int minKSubarraySum(int[] arr, int size){
        int min = 0;
        if (arr.length < size || arr.length == 0){
            return -1;
        }
        //if size is same as array len then sum of array is max/min
        if (arr.length == size){
            for (int i = 0;i < size;i++){
                min += arr[i];
            }
        }
        else{
            int st = 0;
            int end = 1;
            int sum = arr[st];
            //get window size sum first then will traverse array 
            while(end < size){
                sum += arr[end++];
            }
            min = sum;
            //now traverse with window on array 
            while(end < arr.length){
                //changing sum value with incrementing st and end at the same time
                sum -= arr[st++];
                sum += arr[end++];
                if( min > sum){
                    min = sum;
                }

            }

        }

        return min;
    }


    public int minimizedHeight(int[] arr, int k){
        //sort array in nlogn
        Arrays.sort(arr);
        //getting default max managed and min managed value difference
        int result = (arr[arr.length - 1] - k) - (arr[0] + k);
        //checking exceptions in a sorted array(except 1st and last element)
        for(int i = 1; i < arr.length - 1;i++){
            //cause every element has to change(by +/- k)
            int max = Math.max(arr[i] + k, arr[i] - k);
            int min = Math.min(arr[i] + k, arr[i] - k);
            //if found max min difference less change result
            result = Math.min(result, max - min);
        }
        return result;
    }

    public int minJumpsToEnd(int[] arr){
        if (arr.length <= 1) 
            return 0; 
  
        // Return -1 if not possible to jump 
        if (arr[0] == 0) 
            return -1; 
  
        // initialization 
        int maxReach = arr[0]; 
        int step = arr[0]; 
        int jump = 1; 
  
        // Start traversing array 
        for (int i = 1; i < arr.length; i++) { 
            // Check if we have reached  
            // the end of the array 
            if (i == arr.length - 1) 
                return jump; 
  
            // updating maxReach 
            maxReach = Math.max(maxReach, i + arr[i]); 
  
            // we use a step to get to the current index 
            step--; 
  
            // If no further steps left 
            if (step == 0) { 
                // we must have used a jump 
                jump++; 
  
                // Check if the current  
// index/position or lesser index 
                // is the maximum reach point  
// from the previous indexes 
                if (i >= maxReach) 
                    return -1; 
  
                // re-initialize the steps to the amount 
                // of steps to reach maxReach from position i. 
                step = maxReach - i; 
            } 
        } 
  
        return -1; 

    }

    //find duplicate in distinct value array (only 1 value is twice)
    public int duplicate2Val(int[] arr){
        Set<Integer> set = new HashSet<>();
        for (int i = 0;i < arr.length; i++ ){
            //if appraing twice returning that value
            if (set.contains(arr[i])){
                
                return arr[i];
            }
            set.add(arr[i]);
        }
        return -1;
    }
    
    //merge and do inplace chages merge without extra memory so that arr1 < arr2 elements
    public void mergeInplace(int[] arr1, int[] arr2){
        int tail1 = arr1.length - 1;
        int head2 = 0;
        while(tail1 >= 0 && head2 < arr2.length){
            //maintaining smaller values in arr1 via swaping 
            //as 2nd array can have small values in front and
            //arr1 can have bigger values in last as both are sorted
            if(arr1[tail1] > arr2[head2]){
                int temp = arr2[head2];
                arr2[head2] = arr1[tail1];
                arr1[tail1] = temp;
            }
            tail1--;
            head2++;
        }
        //as now arr1 has all the smallest values
        //sorting in nlogn
        Arrays.sort(arr1);
        Arrays.sort(arr2);

    }

    public int[][] mergeIntervals(int[][] intervals){
        LinkedList<int[]> list = new LinkedList<int[]>();
        int maxTime;
        Arrays.sort(intervals, (a, b) -> a[0]-b[0]);
        list.add(intervals[0]);
        for(int i = 1;i < intervals.length;i++){
            if (list.getLast()[1] < intervals[i][0]){
                //if list emplty i.e 1st interval or end time of last interval is less than start time of next intrval
                //add that interval into the list
                list.add(intervals[i]);
            }
            else{
                maxTime = Math.max(list.getLast()[1], intervals[i][1]);
                list.getLast()[1] = maxTime;
            }
        }

        return list.toArray(new int[list.size()][]);
    }
    
    private void reverseNoreturn(int[] arr, int from, int to ){
        int temp;
        //swap from both ends
        while(from < to){
            temp = arr[from];
            arr[from] = arr[to];
            arr[to] = temp;
            from++;
            to--;

        }
    }
    
    public void nextPermutation(int[] arr){
        //find first non increasing number from right
        //find bigger number from that from right 
        //reverse array from first number to last
        int i = arr.length - 2;
        //i till non increasing number from last
        while(i >=0 && arr[i] > arr[i + 1]){i--;}
        //if non increasing number from back is found
        if(i >= 0){
            int j = arr.length - 1;
            //finding 1st number that is greater than that i number from right
            while(j >= 0 && arr[j] < arr[i]){
                j--;
            }
            //swap these numbers
            int temp = arr[j];
            arr[j] = arr[i];
            arr[i] = temp;
        }
        //if not any way revers the array will give sorted array
        reverseNoreturn(arr, i + 1, arr.length - 1);
    }

    private int mergeCount(int[] arr, int[] temp, int l, int mid, int r){
        int i = l;
        int j = mid;
        int k = l;
        int count = 0;
        //i< mid -1 as j starts from mid index na
        while( i <= mid - 1 && j <= r){
            if(arr[i] <= arr[j]){
             temp[k++] = arr[i++];   
            }
            else{
                temp[k++] = arr[j++];
                count += mid - i;
            }
        }
        while(i <= mid - 1){
            temp[k++] = arr[i++];
        }
        while(j <= r){
            temp[k++] = arr[j++];
        }
        //putting back values into arr
        for(i = l;i <= r;i++){
            arr[i] = temp[i];
        }
        //or
        //System.arraycopy(temp, l, arr, l, (r - l + 1));
        show(temp);
        return count;
    }

    private int mergeTweakBreak(int[] arr,int[] temp, int l, int r){
        int count = 0;
        //all index based
        if(l < r){
            int mid = (l+r) >> 1;
            //left array break
            count += mergeTweakBreak(arr, temp, l, mid);
            //right array break
            count += mergeTweakBreak(arr, temp, mid + 1, r);

            //merge time values
            count += mergeCount(arr, temp, l, mid + 1, r);
        }
        return count;
    }

    public int inversionCount(int[] arr){
        //using tweaked merge sort
        int[] temp = new int[arr.length];
        //keeping a temporary array for swaps and putting back the original values to arr
        return mergeTweakBreak(arr, temp, 0, arr.length -1);
    }

    public int bestProfit(int[] prices){
        int lowestBuy = Integer.MAX_VALUE;
        int maxProfit = -1;
        for(int i: prices){
            if(lowestBuy > i){
                lowestBuy = i;
            }
            else{
                if(maxProfit < i - lowestBuy){
                    maxProfit = i - lowestBuy;
                }
            }
        }
        return maxProfit;
    }

    public int countPairSum(int[] arr, int sum){
        int count = 0;
        //using hashmap to store the frequency of numbers
        TreeMap<Integer, Integer> map = new TreeMap<>();
        for(int i = 0; i < arr.length;i++){
            if(map.containsKey(arr[i])){
                map.put(arr[i], map.get(arr[i]) + 1);
            }
            else{
                map.put(arr[i], 1);
            }
        }
        System.out.println(map.keySet());
        for(int i = 0; i < arr.length;i++){
            if(map.containsKey(sum - arr[i])){
                count += map.get(sum - arr[i]);
            }
            //same number twice
            if(sum - arr[i] == arr[i]){
                count--;
            }

        }
        //as we are counting pairs
        return count /2;
    }

    public int maxSubProd(int[] arr){
        int max = 1;
        int min = 1;
        int finalMax = 1;
        for (int i = 0;i < arr.length;i++){
            if(arr[i] == 0){
                finalMax = Math.max(finalMax, 0);    
                max = 1;
                min = 1;
            }
            int a1 = arr[i] * max;
            int a2 = arr[i] * min;
            max = Math.max(Math.max(a1, a2), arr[i]);
            min = Math.min(Math.min(a1, a2), arr[i]);
            finalMax = Math.max(finalMax, max);
        }
        return finalMax;
    }

    public boolean subArrayHasSum0(int[] arr){//O(n)
        //using hasSet
        int sum = 0;
        HashSet<Integer> set = new HashSet<>();
        for(int i = 0;i < arr.length;i++){
            sum += arr[i];
            if(set.contains(sum) || sum == 0 || arr[i] == 0){
                return true;
            }
            else{
                set.add(sum);
            }

        }
        return false;
    } 

    public List<Integer> commonElementIn3Sorted(int[] a,int[] b, int[] c)//T_O(n+m+l)
    {
        List<Integer> list = new ArrayList<Integer>();
        int i = 0;
        int j = 0;
        int k = 0;
        while(i < a.length || j < b.length || k < c.length){
        if(a[i] == b[j] && b[j] == c[k]){
            list.add(a[i]);
            i++;
            j++;
            k++;
        }
        else if(a[i] < b[j]){
            i++;
        }
        else if(b[j] < c[k]){
            j++;
        
        }
        else{
            k++;
        }
        
    }
        return list;
    }

    public void rearrangeAlternate(int[] arr){
        int i = 0;
        int j = 0;
        if(arr.length == 1){
            return;
        }
        while(i < arr.length && arr[i] < 0){
            i++;
        }
        //get 1st + value fron right
        j = i;
        //run to right if found -ve value swap it with i and increment i
        while(j < arr.length){
            if(arr[j] < 0){
                swap(arr, i, j);
                i++;
            }
            j++;
        }
        //O(n)
        i = 1;
        j = 0;
        while(j < arr.length){
            if(arr[j] > 0){
                swap(arr, i, j);
                i += 2;
            }
            j++;
        }
        //O(n+n)

    }

    private static void printR(int[] arr, int size){
	    for(int i = size - 1 ;i >= 0;i--){
	        System.out.print(arr[i]);
	    }
    }

    private static int multiply(int[] result, int rsize, int n){
	    int carry = 0;
	    for(int i = 0;i < rsize;i++){
	        int prod = n*result[i] + carry;
	        result[i] = prod % 10;
	        carry = prod / 10;
	    }
	    while(carry != 0){
	        result[rsize] = carry % 10;
	        carry = carry/10;
	        rsize++;
	        
	    }
	    return rsize;
	}
	
	
	public void longFactorial(int n){
	    int[] result= new int[500];
	    int rsize = 1;
	    result[0] = 1;
	    for(int i = 2;i <= n ;i++){
	        rsize = multiply(result, rsize, i);
	    }
	        printR(result, rsize);
	}
    
    public int lenLongConsecutiveSub(int[] arr){//with out changing array
        //array is unsorted may be rotated
        int len = 0;
        //get 1st index
        int i = 0;
        while(arr[i]  != (arr[i+1] - 1) && i < arr.length -1){
            i++;
        }
        System.out.println(i);
        //found 1st index
        if(i < arr.length -1){
            len++;
            System.out.println(i);
            //it will break as soon as it find non consecutive number
            //if its a sorted array it will not be a infi loop as at end of array last element is compared to +1
            while(arr[i] + 1 == arr[i+1] ){
                len++;
                System.out.println(i);    
                //if goes beyond max index start from 0
                if(i == arr.length-1 && arr[arr.length-1] + 1 == arr[0]){
                    System.out.println(i);    
                    i = 0;
                    len++;
                    
                }else{
                    i++;
                }
                
            }

        }
        return len;
    }

    public static class NodeInt{
        int val;
        int count = 0;
    }
   
    public void freqKMore(int[] arr, int k){//T_O(nK)
        //tetric approch or can implement using simple hashMapin O(n)
        NodeInt[] tetric =  new NodeInt[k];
        //initiateting tectric
        for(int i = 0;i < k;i++){
            tetric[i] = new NodeInt();
        }
        System.out.print("Done initiating..");
        for(int i = 0; i < arr.length;i++){
            System.out.print("In loop iterating..");
        
            int j = 0;
            int isthere = 0;
            //add new element in tetric
            while(j < k){
                if(tetric[j].val == arr[i]){
                    tetric[j].count += 1;
                    System.out.print("found value incrmenting count..");
        
                    break;
                }
                j++;
            }
            //if not in tetric
            if(j == k){
                System.out.print("NOt found in tetric..");
        
                //is there a vacancy
                while(isthere < k){
                    if(tetric[isthere].count == 0){
                        //add new element
                        System.out.print("found 0 count of element.. so adding");
        
                        tetric[isthere].val = arr[i];
                        tetric[isthere].count = 1;
                        break;
                    }
                    isthere++;
                }
            }
            //tetric is full so -=1 to all
            if(isthere == k){
                System.out.print("tetric is full so reducing ...");
        
                isthere = 0;
                while(isthere < k){
                    tetric[isthere].count -=1;
                    isthere++;
                }
                //adding new number if older number ==0, after reducing
                isthere = 0;
                while(isthere < k){
                    if(tetric[isthere].count == 0){
                        tetric[isthere].val = arr[i];
                        tetric[isthere].count = 1;
                        break;
                    }
                    isthere++;
                }
            }
        }
        //tetric values
        for(int i = 0;i < k;i++){
            System.out.println(tetric[i].val);
        }
    }

    //max possible profit can be made for stock prices
    public int maxProfitAtAll(int[] prices){
        int max = 0;
        int prof = 0;
        //adding all positive profits
        //count subsequent profit se if max
        for(int i = 1;i< prices.length;i++){
            prof = prices[i] - prices[i - 1];
            if(prof > 0){
                max += prof;
            }
        }
        return max;
    }
    public boolean ifTripletsum(int[] arr, int sum){
        //if there exist 3 number sum ==target
        //sort the array in T_O(nlogn)
        if(arr.length <= 2){
            return false;
        }

        Arrays.sort(arr);
        int nxt = 0;
        int last = arr.length - 1;
        show(arr);
        for(int i = 0;i < arr.length;i++){
            
            //i stay constant and we move nxt and last in sorted array
            //nxt < last so 3 distinct number then only true
            nxt = i + 1;
            last = arr.length - 1;
            //initializing sumgrid for new i, nxt and last
            int sumgrid = arr[nxt] + arr[last] + arr[i];
            while(nxt < last){
                System.out.println(arr[i]+" "+arr[nxt]+" "+arr[last]+" "+sumgrid);
                if(sumgrid < sum){
                    //take out nxt and increase nxt in sorted array
                    sumgrid -= arr[nxt++];
                    sumgrid += arr[nxt];
                }
                else if(sumgrid > sum){
                    //take out last and reduce last
                    sumgrid -= arr[last--];
                    sumgrid += arr[last];
                }
                else{
                    //sumgrid == sum
                    return true;
                }
            } 
        }
        return false;
    }

    public int trapingRain(int[] arr){
        //idea is get max left and right height, and water traped is min(leftH, rightH) - current index
        //having 2 arrays to sotore hax heights
        int water = 0;
        int[] leftMaxH = new int[arr.length];
        int[] rightMaxH = new int[arr.length];
        
        //init array and adding max if heights from left
        leftMaxH[0] = arr[0];
        for(int i = 1;i < arr.length; i++){
            leftMaxH[i] = Math.max(leftMaxH[i - 1], arr[i]);
        }

        //init array and adding max if heights from right
        rightMaxH[arr.length - 1] = arr[arr.length - 1];
        for(int i = arr.length - 2;i >= 0; i--){
            rightMaxH[i] = Math.max(rightMaxH[i + 1], arr[i]);
        }

        //now traversing the array of heights to find total water, 
        //min(height of left or right at that point) - current height while traversal till that point
        for(int i = 0;i <arr.length; i++){
            water += Math.min(leftMaxH[i], rightMaxH[i])  - arr[i];
        }
        return water;
    }

    public int countMinSubForSumGreaterThanX(int[] arr, int x){//T_O(n)
        //first get the grid sum > x by increasing end then reduce the st and keep trak of
        //min of end - st
        if(arr.length <= 1){
            return -1;
        }
        int min = arr.length;
        int st = 0;
        int end = 1;
        int sum = arr[st];
        if(sum >= x){
            return 1;
        }
        while(end < arr.length){
            while(sum <= x && end < arr.length){
                sum += arr[end++];
            }
        //found optimal end
        //now reducing st
            while (sum > x && st < arr.length){
                if(end - st < min){
                    min = end - st;
                }
                sum -= arr[st++];
                
            }
            
         }
         return min;
    }
    
    //ERROR//////////////////////////////////////////////////////
    /*
    public void partationInto3(int[] arr, int a, int b){
        int last = arr.length - 1;
        int st = 0;
        for(int i = 0;i < last;i++){
            //traversing and swaping less than a values to front
            if(arr[i] <= a){
                int temp = arr[st];
                arr[st] = arr[i];
                arr[i] = temp;
                st++;
            }
            //traversing and swaping greater than b values to last
            else if(arr[i] >= b){
                int temp = arr[last];
                arr[last] = arr[i];
                arr[i] = temp;
                last--;
            }
        }
    }

    public void spiralTraversal(int[][] matrix){
        int rl = matrix.length - 1;
        int cl = matrix[0].length - 1;
        int i = 0;
        int j = 0;
        int curl = 0;
        int s = Math.min(rl+1, cl+1);
        while(curl < s/2 ){
            i = curl;
            j = curl;
            while(j <= cl){
                System.out.println(matrix[i][j++]);
            }
            j--;
            //next row
            if(i + 1 <= rl){i++;}
            while(i <= rl){
                System.out.println(matrix[i++][j]);
            }
            i--;
            if(j - 1 >= curl){j--;}
            while(j >= curl){
                System.out.println(matrix[i][j--]);            
            }
            j++;
            curl++;
            while(i > curl){
            //add matrix[i][j];
            System.out.println(matrix[--i][j]);
            }
            rl--;
            cl--;
        }
        
    }
    */

    
}
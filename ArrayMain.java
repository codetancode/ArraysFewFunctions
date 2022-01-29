import java.util.*;




//Its a new class i am creating inside ArrayMain class, 
//which is implementing Comparable of type new class.
//1.we need to implement compareTo() method for new class
//comapreTo() method is an Abstract method from "Comparable" interface, so new class implements Comparable 
class StrVal implements Comparable< StrVal > {
  int number;//data member 
  //constructor
  public StrVal(int number) {
      this.number = number;
  }
  //@Override- Just notation when we are overriding/implementing abstract method(good practice)
  @Override
  public int compareTo(StrVal o) {
//  	concatinating both number string so that they are of == lengths

      String first = String.valueOf(this.number) + String.valueOf(o.number);
      String second = String.valueOf(o.number) + String.valueOf(this.number);
      //String.compareTo(return difference of 1st char), if same then difference of length of string
//      ex-, "1239".compareTo("91239"))->1-9=(-8)
//      ex-, "91239".compareTo("123"))->9-1=(8)
//      ex-, "1".compareTo("12"))->1-2=(-1)
//      ex-, "1".compareTo("12345"))->1-5=(-4)
      
//  	then comparing with string comparator, based on 1st character
      return second.compareTo(first);
  }
}


public class ArrayMain {
	
	static boolean checkPoss(int diff, int cows, int[] A){
        int putcow=1;
        int recently_placed_cow = A[0];//putting at start
        for(int i=1;i < A.length;i++){
            //check if adjecent differen is there then add cows
            if(A[i] - recently_placed_cow >= diff){
                putcow++;
                //updating as we put a cow
                recently_placed_cow = A[i];
            }
        }
        return putcow >= cows;
    }

    static int maxDiffCow(int[] A, int B) {
    //fixing low and high possibilities range of possible answers
    //lowest difference can be 1 (all index of (linear sorted arrayi.e 1 2 3..) has cows)
    //highesst difference(2 cows at 1st and last index) can be last - 1st
    int low = 1;
    int high = A[A.length-1] - A[0];
    int mid = 0;
    while(low <= high){
        mid = (low + high)>>1;
        if(checkPoss(mid, B, A)){
            //if possible look for better ans go right
            low = mid + 1;
        }
        else{
            //no possible arrangement go left
            high = mid - 1;
        }
    }
    return mid;

    }
    static int sumTill(ArrayList<Integer> a, int till){
        int s=0;
        for(int i=0;i < till;i++){
            s += a.get(i);
        }
        return s;
    }
    
    static boolean checkPainter(ArrayList<Integer> C, int maxTime, int maxGivenPainters, int BtimePerUnit){
        int time = 0;
        int painter = 1;
        int mod = 1000000007;
        //check how many painter req to complete in maxTime
        for(int i=0;i < C.size();i++){
            //assume 1 painter, keep increasing time for additiona boats with given rate (BtimePerUnit)
            time = time + C.get(i)*BtimePerUnit%mod;
            
            //if adding new boat increases the time from maxTime, add a painter to do that new boat
            if(time > maxTime){
                //add new painter
                painter++;
                //now total time will reset to 
                time = C.get(i)*BtimePerUnit%mod;
            }
        }
        //then check if req painter are less than given painter?
        //if req painter are less or == given ->true, else false
        return painter <= maxGivenPainters;

    }
    
    static int solve(String A, String B) {
        HashMap<Character, Integer> toComp = new HashMap<>();
        HashMap<Character, Integer> everyLen = new HashMap<>();
        int tocomplen = A.length();
        int res=0;
        for(int i=0;i < tocomplen;i++){
            if(toComp.containsKey(A.charAt(i))){
                toComp.put(A.charAt(i), toComp.get(A.charAt(i)) + 1);
            }
            else{
                toComp.put(A.charAt(i), 1);
            }
        }
        int i=0, j=0;
        
        while(j < B.length()){
            if(everyLen.containsKey(B.charAt(j))){
                everyLen.put(B.charAt(j), everyLen.get(B.charAt(j)) + 1);
            }
            else{
                everyLen.put(B.charAt(j), 1);
            }
            if( j-i+1 == tocomplen){
                if(toComp.equals(everyLen)){
                    res++;
                }
                else{
                    //remove ith char from everyLen
                    //it will obviously contain char at i, 
                    if(everyLen.containsKey(B.charAt(i))){
                        if(everyLen.get(B.charAt(i)) > 1 ){
                            everyLen.put(B.charAt(i), everyLen.get(B.charAt(i)) - 1);
                        }//if only 1 then remove the key
                        else{
                            everyLen.remove(B.charAt(i));
                        }   
                    } 
                }
                //remove 1 char i and add jth char    
                    i++;
                   
            }

            j++;
        }
        return res;
        }
    
    static String sorted(String tosort){
        int[] charfreq = new int[27];
        String res = "";
        for(int i=0;i < tosort.length();i++){
            charfreq[(int)(tosort.charAt(i)-'a')] ++;
        }
        for(int i=0;i < 27;i++){
            while(charfreq[i] > 0){
                res = res + (char)(i+'a');
                charfreq[i]--;
            }
        }
        return res;
    }
    
    static int gcd(int a, int b){
        if( b == 0)return a;

        return gcd(b, a%b);
    }
    
    static public String minWindow(String A, String B) {
        HashMap<Character, Integer> fix = new HashMap<>();
        HashMap<Character, Integer> checkfix = new HashMap<>();
        
        int[] fromTo = new int[2];
        fromTo[0] = -1;
        fromTo[1] = Integer.MAX_VALUE-2;
        
        for(int i=0;i < B.length();i++){
            fix.putIfAbsent(B.charAt(i), 0);
            fix.replace(B.charAt(i), fix.get(B.charAt(i)) +1);
        }
        int j=0;
        int i=-1;
        while(j < A.length()){
            if(fix.containsKey(A.charAt(j))){
                checkfix.putIfAbsent(A.charAt(j), 0);
                checkfix.replace(A.charAt(j), checkfix.get(A.charAt(j)) +1);
                if(i == -1){
                    i = j;
                }

            }
            if(checkfix.equals(fix)){
                //if diff more update, i, j in ans
                if(fromTo[1] - fromTo[0] > j-i){
                  fromTo[0] = i;
                  fromTo[1] = j;  
                }
                //increase i of window
                if(checkfix.get(A.charAt(i)) > 1){
                    checkfix.replace(A.charAt(i), checkfix.get(A.charAt(i)) - 1);
                }else{
                    checkfix.remove(A.charAt(i));
                }
                i++;
            }
            j++;
        }
        String res =  "";
        if(fromTo[0] != -1 && fromTo[1] != Integer.MAX_VALUE){
            for(int k=fromTo[0];k <= fromTo[1];k++){
                res+=A.charAt(k);
            }
        }
        
        return res;

    }
    static  int solve(String A) {
        int n = A.length();
        int i=0;
        int j=n-1;
        int count =0;
        while(i < j){
            if(A.charAt(i) != A.charAt(j)){
                    count++;
                    j--;
                }else{
                    i++;
                    j--;
                }    
        }
        return count;
    }
    static void permutate(ArrayList<ArrayList<Integer>> store, ArrayList<Integer> arr, int l, int r,int kthPerm){
        if(l == r){
            //new array to store
            ArrayList<Integer> arr1 = new ArrayList<Integer>(arr);
            store.add(arr1);
            if(store.size() == kthPerm) {
            	for(int k =0;k < store.get(kthPerm-1).size();k++) {
            		System.out.print(store.get(kthPerm-1).get(k)+" ");
            	}
            	System.out.println();
            }
            return;
        }
        else{
            //we know we have to swap elements in the given arraylist arr
            for(int i=l; i <= r;i++){
                swap(arr, l, i);
                permutate(store, arr, l+1, r, kthPerm );
                swap(arr, l, i);
            }
        }
       
    }
    static void swap(ArrayList<Integer> arr, int i, int j){
        int t = arr.get(i);
        arr.set(i, arr.get(j));
        arr.set(j, t);
    }
    static ArrayList<ArrayList<Integer>> permute(ArrayList<Integer> A, int kthPerm) {
        ArrayList<ArrayList<Integer>> store = new ArrayList<ArrayList<Integer>>();
        permutate(store, A, 0, A.size()-1, kthPerm);
        return store;
    }
    
    
    static public ArrayList<Integer> solve(ArrayList<ArrayList<Integer>> A) {
        ArrayList<Integer> res = new ArrayList<Integer>();
        Stack<Integer> mainSt = new Stack<>();
        HashMap<Integer, Integer> elementToFreq= new HashMap<>();
        HashMap<Integer, HashSet<Integer>> freqToElement= new HashMap<>();
        int maxFreq = -1;
        //going through query
        for(int i=0;i < A.size();i++){
            if(A.get(i).get(0) == 1){
                // adding to main stack 1st then managing maps
                int nE = A.get(i).get(1);
                //push nE and return -1 so store in res
                mainSt.push(nE); 
                res.add(-1);

                // elementToFreq map
                if(!elementToFreq.containsKey(nE)){
                    //new element will alwasy be freq 1
                    elementToFreq.put(nE, 1);
                    //if very 1st element then, 1 doesnot have a corrosponding set
                    if(!freqToElement.containsKey(1)){
                        freqToElement.put(1, new HashSet<>());
                        freqToElement.get(1).add(nE);
                        //for very first element update maxfreq
                        maxFreq = 1;
                    }
                    else{//just add new element of freq 1
                        freqToElement.get(1).add(nE);
                    }
                }
                else{
                    //update element freq in elementToFreq
                    elementToFreq.replace(nE, elementToFreq.get(nE) + 1);
                    //now go to prev freq set and remove from there and add to new freq set if exist
                    freqToElement.get(elementToFreq.get(nE)-1).remove(nE);
                    int newfreq = elementToFreq.get(nE);
                    //for new freq number add set if not already eist(update max frq), else just add nE to new freq set
                    if(!freqToElement.containsKey(newfreq)){
                        freqToElement.put(newfreq, new HashSet<>());
                        freqToElement.get(newfreq).add(nE);
                        //for very first element update maxfreq
                        maxFreq = newfreq;
                    }
                    else{
                        //just add new element of newfreq
                        freqToElement.get(newfreq).add(nE);
                        if(maxFreq < newfreq){
                            maxFreq = newfreq;
                        }
                    }
                }
            }
            else{//if get maxFreq
                //remove and return max freq element
                //any maxE from max set
                if(maxFreq != -1){

                    Iterator<Integer> itr = freqToElement.get(maxFreq).iterator();
                    int maxE = 0;
                    // if(freqToElement.get(maxFreq).size() > 1){
                    //     while(mainSt.peek() != maxE){
                    //         maxE = itr.next();
                    //     }
                    // }
                    // else{
                        maxE = itr.next();
                    // }
                    res.add(maxE);
                    //now remove maxE from both maps
                    elementToFreq.put(maxE, elementToFreq.get(maxE) - 1);
                    //remove maxE from next frq and add to previous freq
                    freqToElement.get(elementToFreq.get(maxE)+1).remove(maxE);
                    //not removing any set, so set will be present for all freq
                    //for the case of 0 set wont exist so dont add it anything, add only if prev frq exist
                    if(freqToElement.containsKey(elementToFreq.get(maxE))){
                        freqToElement.get(elementToFreq.get(maxE)).add(maxE);
                        //if removed from max freq set
                        if(freqToElement.get(elementToFreq.get(maxE)+1).isEmpty()){
                            maxFreq = maxFreq - 1;
                        }
                    }else{
                        //if element is remove from 1 freq set
                        if(freqToElement.get(1).isEmpty()){
                            freqToElement.remove(1);
                            maxFreq = -1;
                        }
                    }
                    if(elementToFreq.get(maxE) == 0){
                        //remove the key
                        elementToFreq.remove(maxE);
                    }
                    //pop from mainst, 
                    mainSt.pop();  
                }
            }
        }
        return res;
    }
    static public ArrayList<Integer> solveStack(ArrayList<ArrayList<Integer>> A) {
        ArrayList<Integer> res = new ArrayList<Integer>();
        Stack<Integer> mainSt = new Stack<>();
        HashMap<Integer, Integer> elementToFreq= new HashMap<>();
        HashMap<Integer, Stack<Integer>> freqToElement= new HashMap<>();
        int maxFreq = -1;
        //going through query
        for(int i=0;i < A.size();i++){
            if(A.get(i).get(0) == 1){
                // adding to main stack 1st then managing maps
                int nE = A.get(i).get(1);
                //push nE and return -1 so store in res
                mainSt.push(nE); 
                res.add(-1);

                // elementToFreq map
                if(!elementToFreq.containsKey(nE)){
                    //new element will alwasy be freq 1
                    elementToFreq.put(nE, 1);
                    //if very 1st element then, 1 doesnot have a corrosponding set
                    if(!freqToElement.containsKey(1)){
                        freqToElement.put(1, new Stack<>());
                        freqToElement.get(1).push(nE);
                        //for very first element update maxfreq
                        maxFreq = 1;
                    }
                    else{//just add new element of freq 1
                        freqToElement.get(1).push(nE);
                    }
                }
                else{
                    //update element freq in elementToFreq
                    elementToFreq.replace(nE, elementToFreq.get(nE) + 1);
                    int newfreq = elementToFreq.get(nE);
                    
                    if(!freqToElement.containsKey(newfreq)){
                        freqToElement.put(newfreq, new Stack<>());
                        freqToElement.get(newfreq).push(nE);
                        //for very first element update maxfreq
                        maxFreq = newfreq;
                    }
                    else{
                        //just add new element of newfreq
                        freqToElement.get(newfreq).push(nE);
                        if(maxFreq < newfreq){
                            maxFreq = newfreq;
                        }
                    }
                }
            }
            else{//if get maxFreq
                //remove and return max freq element
                //any maxE from max set
                if(maxFreq != -1){
                    int maxE = freqToElement.get(maxFreq).pop();
                   
                    res.add(maxE);
                    //now remove maxE from both maps
                    elementToFreq.put(maxE, elementToFreq.get(maxE) - 1);
                    if(freqToElement.get(elementToFreq.get(maxE)+1).isEmpty()){
                        maxFreq = maxFreq - 1;
                    }
                    //if element is remove from 1 freq set
                    if(freqToElement.get(1).isEmpty()){
                        freqToElement.remove(1);
                        maxFreq = -1;
                    }
              
                    if(elementToFreq.get(maxE) == 0){
                        //remove the key
                        elementToFreq.remove(maxE);
                    }
                    //pop from mainst, 
                    mainSt.pop();  
                }
            }
        }
        return res;
 }
    
    //idea is to find fnal sign for everry expression (a, b, c) as per contraints
    //so there will be local sign and global sign true as positive, false as -
    static boolean localSign(String from, int atI){
        if(atI-1 < 0){
            return true;
        }
        if(from.charAt(atI-1) != '-'){
            return true;
        }
        else{
            return false;
            }
    }
    //for global sign we are maintaining Stack
    static boolean ifAlpha(char test){
        if(test !='-' && test !='+' && test !=')' && test !='(' ){
            return true;
        }
        else{return false;}
    }
    static int solveMatch(String A, String B) {
        HashMap<Character, Boolean> Aplusminus = new HashMap<>();
        HashMap<Character, Boolean> Bplusminus = new HashMap<>();
        Stack<Boolean> globalSign = new Stack<>();
        globalSign.push(true);
        //+ initially
        //if size is > 1 and peek is (, then gbl sut be minus
        for(int i=0;i < A.length();i++){
            char temp = A.charAt(i);
            if(ifAlpha(temp)){
                //update map and get final signs, from global and local
                //update its signs
                if(globalSign.peek()){
                    Aplusminus.put(temp, localSign(A, i));
                }else{
                    Aplusminus.put(temp, !localSign(A, i));
                }
            }
            else{
                //update stack
                if(temp == '('){
                    //update global sign, if + repeate top, if -, not of top
                    if(localSign(A, i)){
                        globalSign.push(globalSign.peek());
                    }else{
                        globalSign.push(!globalSign.peek());
                    }
                    
                }else{//closing breacket so pop 1 sign
                    if(temp == ')'){
                        globalSign.pop();
                    }
                }
            }
        }
        //rset globalSign
        globalSign.clear();
        globalSign.push(true);
        for(int i=0;i < B.length();i++){
            char temp = B.charAt(i);
            if(ifAlpha(temp)){
                //update map and get final signs, from global and local
                //update its signs
                if(globalSign.peek()){
                    Aplusminus.put(temp, localSign(B, i));
                }else{
                    Aplusminus.put(temp, !localSign(B, i));
                }
                
            }
            else{
                //update stack
                if(temp == '('){
                    //update global sign, if + repeate top, if -, not of top
                    if(localSign(B, i)){
                        globalSign.push(globalSign.peek());
                    }else{
                        globalSign.push(!globalSign.peek());
                    }
                    
                }else{//closing breacket so pop 1 sign
                    if(temp == ')'){
                        globalSign.pop();
                    }
                }
            }
        }

        //compare both maps
        if(Aplusminus.equals(Bplusminus)){
            return 1;
        }
        else{
            return 0;
        }
    }
    static ArrayList<Integer> z_fun(String A){
        ArrayList<Integer> zOfA = new ArrayList<Integer>(Collections.nCopies(A.length(), 0));
        int l=0;
        int r=0;
        //l ,r, for window, curr starting from 1 as 0 will be len of string only
        for(int curr = 1; curr < A.length();curr++){
            if(curr > r){
                //bruteforce
            	 l=r=curr;
                //while stays within len, and r-l is same as r(first(r-l)index is same as new r, r++)
                while(r < A.length() && A.charAt(r-l) == A.charAt(r)){r++;}
                //now set curr, to the extent r went(r-l)
                zOfA.set(curr, r-l);
                //putting r inside the window of r and l, by r-- after while loop
                r--;
            }else{
                //curr inside the window of l and r so copy
                //if curr index value + curr-l index value stays <= r copy usual
                if(curr + zOfA.get(curr-l) <= r){
                    zOfA.set(curr, zOfA.get(curr-l));//and curr will change with for loop
                }else{
                    //case where window value exceeds the window size of r
                    //setting curr and l as same and check for new value for curr via bruteforce
                    l=curr;
                    while(r < A.length() && A.charAt(r-l) == A.charAt(r)){r++;}
                    zOfA.set(curr, r-l);
                     //putting r inside the window of r and l, by r-- after while loop
                    r--;
                } 
            }

        }
        return zOfA;
    }
  
    static void showAL(ArrayList<Integer> al) {
    	System.out.print("{ ");
    	for(int i=0;i < al.size();i++) {
    		System.out.print(al.get(i)+" ");
    	}
    	System.out.println(" }");
    }
    static String smallestPrefix(String A, String B) {
        String res = "";
        //lexicographically smallest string of both prexies of Ab and, so both 1st of A and B will be in the answer string
        int j=0;
        if(A.charAt(0) <= B.charAt(0)){
            while(j < A.length() && A.charAt(j) <= B.charAt(0)){
            res+=A.charAt(j);
            j++;
            }
        }else{
            res+=A.charAt(j);
            j++;//|| A.charAt(j) <= B.charAt(0)
            while(j < A.length() && (A.charAt(j) <= A.charAt(0) )){
                res+=A.charAt(j);
                j++;
                
            }
        }
        res+=B.charAt(0);
        

        return res;
    }
    static String revStr(String A){
        int i=0;
        String rev = "";
        while(i < A.length()){
            char temp = A.charAt(i);
            rev=temp+rev;
            i++;
        }
        return rev;
    }
    //z function check
    static ArrayList<Integer> kmp_fun(String A){
        ArrayList<Integer> lps = new ArrayList<Integer>(Collections.nCopies(A.length(), 0));
        //longest prefix suffix
        int i=0;
        int j=1;
        while(j < A.length()){
            if(A.charAt(i) == A.charAt(j)){
                i++;
                //update lps array
                lps.set(j, i);
                j++;
            }else{
                //if a mismatch
                //if i is not at 0 index
                if(i != 0){
                    //i takes lps of i-1 value
                    i=lps.get(i-1);
                }else{
                    //if i at 0 index,update lps of j
                    lps.set(j, 0);
                    //  move j if there is a mismatch
                    j++;
                }
            }

        }
        return lps;
    }
    static int solveStringPal(String A) {
        int n = A.length();
        //append same string with an distinct char in middle
        String concat = A + "$" + revStr(A);
        //now find lps of concat, 
        // difference of stringlen(A) - last index of lps will say how many different char are there.
        ArrayList<Integer> lps = kmp_fun(concat);
        showAL(lps);
        return (lps.get(lps.size()-1)-n);

        }
    static int evalRPN(ArrayList<String> A) {
        int resultant = 0;
        Stack<String> st = new Stack<>();
        HashSet<String> operator = new HashSet<>();
        operator.add("*");
        operator.add("+");
        operator.add("-");
        operator.add("/");
        
        int i=0;
        while(i < A.size()){
            String temp = A.get(i);
            if(!st.isEmpty() && operator.contains(temp)){
                int op2 = Integer.parseInt(st.pop());
                int op1 = Integer.parseInt(st.pop());
                if(temp == "*"){
                    resultant = op1*op2;
                }
                else if(temp.equals("/") ){
                    resultant = op1/op2;
                }
                else if(temp.equals("+") ){
                    resultant = op1+op2;
                }
                else if(temp.equals("-") ){
                    resultant = op1-op2;
                }
                //finally push
                st.push(Integer.toString(resultant));
            }else{
                // if not operator push numbers
                st.push(temp);
            }
            i++;
        }

        return Integer.parseInt(st.pop());
    }
    static  ArrayList<Integer> nextGreater(ArrayList<Integer> A) {
        //just so we dont have to reverse res.
        ArrayList<Integer> res = new ArrayList<Integer>(Collections.nCopies(A.size(), 0));
        Stack<Integer> st = new Stack<>();
        //maintain a stack storing large vaalues goinf from right to left
        int n = A.size()-1;
        for(int i=n;i >= 0;i--){
            //if not not empty and top element is less/= to new element from right to left pop 
            if(!st.isEmpty() && st.peek() <= A.get(i)){
                st.pop();
            }
            //after pop, if emty or the right greater element 
            if(st.isEmpty()){
                res.set(i, -1);
            }else{
                res.set(i, st.peek());
            }
            //anyway push new element for futute comparisions
            st.push(A.get(i));
        }
        return res;
    }
    static public int solveMAXMIN(ArrayList<Integer> A) {
        ArrayList<Integer> smallerIndexFromLeft = new ArrayList<Integer>();
        ArrayList<Integer> smallerIndexFromRight = new ArrayList<Integer>();
        ArrayList<Integer> greaterIndexFromLeft = new ArrayList<Integer>();
        ArrayList<Integer> greaterIndexFromRight = new ArrayList<Integer>();
        Stack<Integer> st = new Stack<>();
        long maxSum = 0;
        long minSum = 0;
        
        int mod=1000000007;
        int n = A.size()-1;
        // filling 2 smaller
        //first left to right, storing index of smaller value found
        for(int i=0;i <= n;i++){
            //if not empty st, and st top has index of greater element, pop()
            while(!st.isEmpty() && (long)A.get(st.peek()) >= (long)A.get(i) ){
                st.pop();
            }
            if(st.isEmpty()){
                smallerIndexFromLeft.add(-1);
            }else{
                smallerIndexFromLeft.add(st.peek());
            }
            //storing index
            st.push(i);
        }

        //clear stack()
        st.clear();

        //second right to left, storing index of smaller value found
        for(int i=n;i >= 0;i--){
            //if not empty st, and st top has index of greater element, pop()
            while(!st.isEmpty() && (long)A.get(st.peek()) >= (long)A.get(i) ){
                st.pop();
            }
            if(st.isEmpty()){
                smallerIndexFromRight.add(0, -1);
            }else{
                smallerIndexFromRight.add(0, st.peek());
            }
            //storing index
            st.push(i);
        }
        
        //clear stack()
        st.clear();
        
        // filling 2 greater 
        //first left to right, storing index of greater value found
        for(int i=0;i <= n;i++){
            //if not empty st, and st top has index of smaller element, pop()
            while(!st.isEmpty() && (long)A.get(st.peek()) <= (long)A.get(i) ){
                st.pop();
            }
            if(st.isEmpty()){
                greaterIndexFromLeft.add(-1);
            }else{
                greaterIndexFromLeft.add(st.peek());
            }
            //storing index
            st.push(i);
        }

        //clear stack()
        st.clear();

        //second right to left, storing index of greater value found
        for(int i=n;i >= 0;i--){
            //if not empty st, and st top has index of smaller element, pop()
            while(!st.isEmpty() && (long)A.get(st.peek()) <= (long)A.get(i) ){
                st.pop();
            }
            if(st.isEmpty()){
                greaterIndexFromRight.add(0, -1);
            }else{
                greaterIndexFromRight.add(0, st.peek());
            }
            //storing index
            st.push(i);
        }


        // finally, contribution technique, 
        // (no. of sub array in which A[i] is max - no. of sub array in which A[i] is max)*A[i]-maxsum
        // (no. of sub array in which A[i] is min - no. of sub array in which A[i] is min)*A[i]-minsum
        for(int i=0;i <= n;i++){
            long sOnLeft = smallerIndexFromLeft.get(i) != -1 ?smallerIndexFromLeft.get(i) :0;
            long sOnRight = smallerIndexFromRight.get(i) != -1 ?smallerIndexFromRight.get(i) :0;
            long gOnLeft = greaterIndexFromLeft.get(i) != -1 ?greaterIndexFromLeft.get(i) :0;
            long gOnRight = greaterIndexFromRight.get(i) != -1 ?greaterIndexFromRight.get(i) :0;

            //(no. of subarr, in which element is max or min) * value
            maxSum = maxSum%mod + ( ( ((i-sOnLeft) * (sOnRight-i) )%mod)*A.get(i) )%mod; 
            minSum = minSum%mod + ( ( ((i-gOnLeft) * (gOnRight-i) )%mod)*A.get(i) )%mod;
        }
        return (int)((maxSum - minSum)%mod);
    }
    static int solveQ(ArrayList<Integer> A, int B) {
        //similar to sliding window, instead of window, use q for add and del elemenents
        Queue<Integer> window= new LinkedList<>();
        long max = Integer.MIN_VALUE;
        long min = Integer.MAX_VALUE;
        long sumMax=0;
        long sumMin=0;
        int mod = 1000000007;
        //first add b elements then drag window through the ArrayList
        int i=0;
        while(i < B){
            max = Math.max(max, A.get(i));
            min = Math.min(min, A.get(i));
            window.add(A.get(i));
            i++;
        }
        sumMax = max;
        sumMin = min;
        while(i < A.size()){
            //reduce sum, remove from window then add sum , add in window then update answer
            int windowout = window.poll(); 
            if(windowout == max){
                //reset max
                max = Integer.MIN_VALUE;
            }
            if(windowout == min){
                //reset min
                min = Integer.MAX_VALUE;
            }
            //add new elements, update with new element
            //updataing ans
            max = Math.max(max, A.get(i));
            min = Math.min(min, A.get(i));
            
            //update sum as new grid of b elements is formed
            sumMax= sumMin%mod + max%mod;
            sumMin= sumMin%mod + min%mod;
            //now add new element in window
            window.add(A.get(i));
            
            i++;
        }


        return (int)(sumMax%mod + sumMin%mod)%mod;
    }

    	
    	
    static ArrayList<Integer> solvePosNeg(ArrayList<Integer> A) {
        ArrayList<Integer> res = new  ArrayList<Integer>();
        int zero=0;
        int one=0;
        while(zero < A.size() && one < A.size()){
            while(zero < A.size() && A.get(zero) >= 0 ){
                zero++;
            }
            if(zero < A.size()){
                res.add(A.get(zero));
                zero++;
            }else{
                break;
            }
            while(one < A.size() && A.get(one) < 0 ){
                one++;
            }
            if(one < A.size()){
                res.add(A.get(one));
                one++;
            }else{
                break;
            }
        }
        
        while(one < A.size() && res.size() != A.size()){
        	System.out.println(zero+" pos->"+A.get(zero));
                
            res.add(A.get(one));
            one++;
        }
        while(zero < A.size() && res.size() != A.size()){
        	System.out.println(zero+" neg->"+A.get(zero));
            )
            res.add(A.get(zero));
            zero++;
        }
        return res;
    }
    
    static void largestRectangleArea(ArrayList<Integer> A) {
        Stack<Integer> st = new Stack<>();
        ArrayList<Integer> leftMin = new ArrayList<Integer>(Collections.nCopies(A.size(), 0));
        ArrayList<Integer> rightMin = new ArrayList<Integer>(Collections.nCopies(A.size(), 0));
        int n = A.size();
        for(int i=0;i < n;i++){
            while(!st.isEmpty() && A.get(st.peek()) > A.get(i)){
                st.pop();
            }
            if(st.isEmpty()){
                leftMin.set(i, -1);
            }
            else{
                leftMin.set(i, st.peek());
            }
            st.push(i);
        }
        st.clear();
        for(int i=n-1;i >=0 ;i--){
            while(!st.isEmpty() && A.get(st.peek()) > A.get(i)){
                st.pop();
            }
            if(st.isEmpty()){
                rightMin.set(i, -1);
            }
            else{
                rightMin.set(i, st.peek());
            }
            st.push(i);
        }
    
        showAL(leftMin);
        showAL(rightMin);
    }
    
    
    
    static int solveMaxMin(ArrayList<Integer> A) {
        ArrayList<Integer> smallerIndexFromLeft = new ArrayList<Integer>(Collections.nCopies(A.size(),0));
        ArrayList<Integer> smallerIndexFromRight = new ArrayList<Integer>(Collections.nCopies(A.size(),0));
        ArrayList<Integer> greaterIndexFromLeft = new ArrayList<Integer>(Collections.nCopies(A.size(),0));
        ArrayList<Integer> greaterIndexFromRight = new ArrayList<Integer>(Collections.nCopies(A.size(),0));
        Stack<Integer> st = new Stack<>();
        long sum = 0;
        
        int mod=1000000007;
        int n = A.size()-1;
        // filling 2 smaller
        //first left to right, storing index of smaller value found
        for(int i=0;i <= n;i++){
            //if not empty st, and st top has index of greater element, pop()
            while(!st.isEmpty() && A.get(st.peek()) > A.get(i) ){
                st.pop();
            }
            if(st.isEmpty()){
                smallerIndexFromLeft.set(i, -1);
            }else{
                smallerIndexFromLeft.set(i, st.peek());
            }
            //storing index
            st.push(i);
        }

        //clear stack()
        st.clear();

        //second right to left, storing index of smaller value found
        for(int i=n;i >= 0;i--){
            //if not empty st, and st top has index of greater element, pop()
            while(!st.isEmpty() && (long)A.get(st.peek()) >= (long)A.get(i) ){
                st.pop();
            }
            if(st.isEmpty()){
                smallerIndexFromRight.set(i, -1);
            }else{
                smallerIndexFromRight.set(i, st.peek());
            }
            //storing index
            st.push(i);
        }
        
        //clear stack()
        st.clear();
        
        // filling 2 greater 
        //first left to right, storing index of greater value found
        for(int i=0;i <= n;i++){
            //if not empty st, and st top has index of smaller element, pop()
            while(!st.isEmpty() && (long)A.get(st.peek()) < (long)A.get(i) ){
                st.pop();
            }
            if(st.isEmpty()){
                greaterIndexFromLeft.set(i, -1);
            }else{
                greaterIndexFromLeft.set(i, st.peek());
            }
            //storing index
            st.push(i);
        }

        //clear stack()
        st.clear();

        //second right to left, storing index of greater value found
        for(int i=n;i >= 0;i--){
            //if not empty st, and st top has index of smaller element, pop()
            while(!st.isEmpty() && (long)A.get(st.peek()) <= (long)A.get(i) ){
                
                st.pop();
            }
            if(st.isEmpty()){
                greaterIndexFromRight.set(i, -1);
            }else{
                greaterIndexFromRight.set(i, st.peek());
            }
            //storing index
            st.push(i);
        }
        showAL(smallerIndexFromLeft);
        showAL(smallerIndexFromRight);
        showAL(greaterIndexFromLeft);
        showAL(greaterIndexFromRight);
        

        // finally, contribution technique, 
        // (no. of sub array in which A[i] is max - no. of sub array in which A[i] is max)*A[i]-maxsum
        // (no. of sub array in which A[i] is min - no. of sub array in which A[i] is min)*A[i]-minsum
        System.out.println("For array this:");
        showAL(A);
        for(int i=0;i <= n;i++){
            long sOnLeft = i - smallerIndexFromLeft.get(i);
            long sOnRight = (smallerIndexFromRight.get(i) != -1 ?smallerIndexFromRight.get(i) :(n+1)) - i;
            long gOnLeft = i - greaterIndexFromLeft.get(i);
            long gOnRight = (greaterIndexFromRight.get(i) != -1 ?greaterIndexFromRight.get(i) :(n+1)) - i;
            
            System.out.println("Sub arrays of "+A.get(i)+" as minimum are  "+(sOnLeft)+" on left.");
            System.out.println("Sub arrays of "+A.get(i)+" as minimum are "+(sOnRight)+" on right.");
            
            System.out.println("Sub arrays of "+A.get(i)+" as maximum are  "+(gOnLeft)+" on left.");
            System.out.println("Sub arrays of "+A.get(i)+" as maximum are "+(gOnRight)+" on right.");
            
            //(no. of subarr, in which element is max or min) * value
            System.out.println("Contribution of "+A.get(i)+" as minimum in "+(sOnLeft) * (sOnRight)+"sub arrays");
            System.out.println("Contribution of "+A.get(i)+" as maximun in "+(gOnLeft) * (gOnRight)+"sub arrays");
            
            sum = sum%mod + ( ( ((gOnLeft) * (gOnRight) )%mod)*A.get(i) )%mod;
            sum = sum%mod - ( ( ((sOnLeft) * (sOnRight) )%mod)*A.get(i) )%mod; 
            if(sum < 0) {
            	sum+=mod;
            }
            sum=sum%mod;
            
        }
        return (int)(sum);
    }
    static void test(ArrayList<Integer> t) {
    	int i=0;
    	int n=t.size();
    	int j=n-1;
    	while(i < n - 1 && t.get(i) < t.get(i+1) ){
    		i++;
    	}
    	while(j >= 0  && t.get(j-1) < t.get(j) ){
    		j--;
    	}
    	System.out.println(i);
    	System.out.println(j);
    	int min = 1000;
    	int max = 1;
    	for(int k=i;k <= j;k++) {
    		min= Math.min(min, t.get(k));
    		max= Math.max(max, t.get(k));
    	}
    	int l=0;
    	int r=n-1;			
    	while(l <= i && t.get(l) <= min  ) {
    		l++;
    	}					
    	while(r >= j && t.get(r) >= max ) {
    		r--;
    	}
    	System.out.println(l);
    	System.out.println(r);
    	
    }
    
    static int solveMax(ArrayList<Integer> A, int B) {
        int minval = Integer.MAX_VALUE;
        int maxval = Integer.MIN_VALUE;
    
        for(int i=0;i < A.size();i++){
            minval = Math.min(minval, A.get(i));
            maxval = Math.max(maxval, A.get(i));
        }

        ArrayList<Integer> freq = new ArrayList<Integer>(Collections.nCopies(maxval+1, 0));
        for(int i=0;i < A.size();i++){
            freq.set(A.get(i), freq.get(A.get(i)) + 1);
        }

        while(B > 0 && maxval > minval ){
            // min++
            if(freq.get(minval) < freq.get(maxval)){
                if(freq.get(minval) < B){
                    int lastMinfreq = freq.get(minval);
                    B = B - freq.get(minval);
                    minval++;
                    freq.set(minval, freq.get(minval) + lastMinfreq);
                }else{
                    break;
                }
            }else{
                // max--
                if(freq.get(maxval) < B){
                    int lastMaxfreq = freq.get(maxval);
                    B = B - freq.get(maxval);
                    maxval--;
                    freq.set(maxval, freq.get(maxval) + lastMaxfreq);
                }else{
                    break;
                }
            }

        }
        
        return (int)(maxval-minval);
    }
    
    
    static int solveMINJ(int[] A) {
        //edge case, no jumps possible
       if(A[0]==0){return -1;}
//  3 2 4 5 5
        int steps = A[0];
        int maxreach = A[0];
        int jumps = 1;
        int n = A.length;
        // [ 1  ]
        
        for(int i=1;i < n;i++){
            if(i == n-1){
                return jumps;
            }                            
            maxreach = Math.max(maxreach, i + A[i]);//1+2(3),3+4(7)5+4(9) 
            //reducing steps
            steps--;//0103

            if(steps == 0){
                jumps++;//2(3)  
                //if we -ve numbers or low numbers 0 then maxreach == 1
                if(i == maxreach){
                    return -1;
                }
                steps = maxreach - i;//3-1 (2)(7-3)(4)
            }
        }
        return -1;
    }
    
    static String solveNonR(String A) {
        Queue<Character> charQ = new LinkedList<>();
        HashMap<Character, Integer> repeatedChar = new HashMap<>();
        StringBuilder topAns = new StringBuilder();

        for(int i=0;i < A.length();i++){
            if(!charQ.isEmpty()){
                char newChar = A.charAt(i);
                
                if(charQ.peek() == newChar ){
                    charQ.poll();   
                    repeatedChar.put(newChar, repeatedChar.get(newChar) + 1);
                }else{
                    //increasing freq of new emelents and adding it
                	if(!repeatedChar.containsKey(newChar)){
                        repeatedChar.put(newChar, 1);
                        charQ.add(newChar);
                    }else {
                    	repeatedChar.put(newChar, repeatedChar.get(newChar) + 1);
                    }
                    
                }

                //now update answer (topans)
                if(!charQ.isEmpty()){
                    while(!charQ.isEmpty() && repeatedChar.get(charQ.peek()) > 1){
                        charQ.poll();
                    }
                    System.out.print("After polling >running  char"+i);
                    System.out.println(charQ.peek());
                    
                    if(!charQ.isEmpty() && repeatedChar.get(charQ.peek()) == 1){
                        topAns.append(charQ.peek());
                    }else {
                        topAns.append("#");
                    }    
                }else{
                    topAns.append("#");
                }    
            }else{
                //update map
            	 if(!repeatedChar.containsKey(A.charAt(i))){
                     repeatedChar.put(A.charAt(i), 1);
                 }
            	 else {
            	     repeatedChar.put(A.charAt(i), repeatedChar.get(A.charAt(i)) + 1);	 
            	 }
            	 //update answer and q
            	charQ.add(A.charAt(i));
            	if(topAns.isEmpty()) {
            		topAns.append(A.charAt(i));
            	}else {
            		topAns.append("#");
            	}
                
            }
        }
        showQueue(charQ);
        return topAns.toString(); 
    
    }
    
    static void showQueue(Queue<Character> a) {
    	System.out.println();
    	System.out.print("[ ");
    	while(a.size() > 0) {
    		System.out.print(a.poll()+" ");
    	}
    	System.out.print(" ]");
    	System.out.println();
    }
    
    static int conSum(int A) {
        ArrayList<Integer> conNum = new ArrayList<>();
        int sum = 0;
        int count = 1;

        for(int temp = 1;sum < A;temp++){
            conNum.add(temp);
            sum+=temp;
        }
        int i = 0;
        int j = conNum.size()-1;
        // [1 , 2]
        System.out.println("first arraylist ready");
        while(i < j && conNum.get(i) < A){
        	
            if(sum == A){
                count++;
                conNum.add(conNum.get(j)+1);
                j++;
                sum += conNum.get(j);
            }
            if(sum > A){
                sum -= conNum.get(i);
                i++;
            }else{
                conNum.add(conNum.get(j)+1);
                j++;
                sum += conNum.get(j);
            }
            
        }
        showAL(conNum);
        return count;
  }
    
    
    
    
    
    static int res = 0;
    static void permutate(ArrayList<Integer> given, int currIndex){
        if(currIndex == given.size()-1){
            if(ifPerfSqare(given.get(currIndex) + given.get(currIndex-1))){
                res++;
                showAL(given);
            }
            return;
        }
        else{
            //we know we have to swap elements in the given arraylist arr
            //keeping swaping and unswap(back tracking) from left to right
            for(int i=currIndex; i <= given.size()-1;i++){
                if(given.get(currIndex) != given.get(i) || i == currIndex){
                    swap(given, currIndex, i);
                    if(currIndex == 0 || ifPerfSqare(given.get(currIndex-1) + given.get(currIndex))){
                            permutate(given, currIndex + 1);
                    }
                    swap(given, currIndex, i);
                }
            }
        }
    }
    
    static boolean ifPerfSqare(int a){
        int res1 = (int)Math.sqrt(a);
        return res1*res1 == a;
    }
    
    
    
    static double findMedianSortedArrays(final List<Integer> a, final List<Integer> b) {//keeping a min
    if(a.size() > b.size()){ 
        return findMedianSortedArrays(b, a);
    }
    //for edge case
    if(a.size() == 0){
        int onlyArrlen = b.size();
        if(onlyArrlen == 1){
            return (double)b.get(0);
        }else{
            if( onlyArrlen%2 == 0){
                // System.out.println( ( b.get( (onlyArrlen/2) ) + b.get( ((onlyArrlen/2)-1) ) )/2  );
                return (double)( ( b.get( (onlyArrlen/2) ) + b.get( ((onlyArrlen/2)-1) ) )/2 );
            }else{
                return (double)(b.get(onlyArrlen/2));
            }
        }
    }
    
    //a always refer to small size array and we put 2 pointer on A, to partation for overall left of mergeArray
    int leftOfAstart = 0;
    int leftOfAend = a.size();
    int midOfmergedArraylen = (a.size() + b.size()+1)/2;
    //so partation from a and b should sum up to mid midOfmergedArraylen
    while(leftOfAstart <= leftOfAend){
        int mid = (leftOfAstart + leftOfAend)>>1;
        int ApartSize = mid;
        int BpartSize = midOfmergedArraylen - mid;

        // if a and b sizes are correct, then partationare just size and size-1 index values
        // and partation in A is Aleft and Aright, similarly partation in B is bleft and b rit
        //other wise min value(just so we dont update answer)
        // part from A
        int leftA = (ApartSize > 0 )?a.get(ApartSize-1):Integer.MIN_VALUE;
        int rightA = (ApartSize < a.size() )?a.get(ApartSize):Integer.MAX_VALUE;
        // part from B
        int leftB = (BpartSize > 0 )?b.get(BpartSize-1):Integer.MIN_VALUE;
        int rightB = (BpartSize < b.size() )?b.get(BpartSize):Integer.MAX_VALUE;
        
        // if left of A value<= right of B and left of b <= right of A, that means both 4 values will be a part of left of merged array
        // so then correct paration for median
        if(leftA <= rightB && leftB <= rightA){
            //now if the merged length is even then mediean is rightA+leftB/2
            // else max of left a or leftB
            if((a.size()+b.size())%2 == 0){
                double maxfromA = Math.max(leftA, rightA);
                double minfromB = Math.min(leftB, rightB);
                System.out.println(leftA+"  "+leftB+", "+rightA+"  "+rightB);
                return (double)((maxfromA+minfromB)/2);
                
            }else{
                System.out.println(leftA+"  "+leftB+", "+rightA+"  "+rightB);
                return (double)Math.max(leftA, leftB);
            }
        }
        else if(leftA > rightB){
            leftOfAend = mid - 1;
        }else{
            leftOfAstart = mid + 1;
        }
    }
    return 0;

	}
    
    
    static int smallerNumThan(int target, ArrayList<ArrayList<Integer>> A){
        int count = 0;
        for(int i=0;i < A.size();i++){
            int low = 0;
            int high = A.get(0).size()-1;
            int mid = 0;
            int temp=0;
            while(low <= high){
                mid = (low+high)>>1;
                if(A.get(i).get(mid) <= target){
                    low=mid+1;
                    temp = mid+1;
                }else{
                    high=mid-1;
                }
            }
            //as 0 based index
            count+=(temp);
        }
        // System.out.println(count);
        return count;
    }
    static int findMedian(ArrayList<ArrayList<Integer>> A) {
        //idea is to find min and max of from 1st and last col of matrix and do BS on it-O(row)
        //checkfunction - check the count of smaller numbers in every row
        //we need to find (1+(n*m)/2)th smaller value(if present that's the median)
        int low = 0;
        int high = Integer.MAX_VALUE;
        int mid =0;
        int desiredCount = 0;
        int ans = 1;
        desiredCount = (  1 + ( (A.size()*A.get(0).size()))/2) ;
        while(low < high){
            mid = (low+high)>>1;
            if(smallerNumThan(mid, A) < desiredCount ){
                low = mid+1;
                
            }else{
                 ans = mid;
                high = mid;
            }
        }
        return ans;

    }
    
    
    static int count = 0;
    static  int numSquarefulPerms(ArrayList<Integer> A) {
        if (A.size() == 0) {
            return 0;
        }
        
        Collections.sort(A);
        boolean[] used = new boolean[A.size()];
        backtracking(A, used, new ArrayList<>());
        return count;
    }
    static void backtracking(ArrayList<Integer> A, boolean[] used, List<Integer> list) {
        if (list.size() == A.size()) {
            count ++;
            return;
        }
        for (int i = 0; i < A.size(); i ++) {
            if (used[i] || (i - 1 >= 0 && A.get(i) == A.get(i-1) && !used[i - 1])) continue;
            if (list.size() == 0 || isSquare(list.get(list.size() - 1) + A.get(i))) {
                list.add(A.get(i));
                used[i] = true;
                backtracking(A, used, list);
                list.remove(list.size() - 1);
                used[i] = false;
            }
        }
    }
    static boolean isSquare(int x) {
        double r = Math.sqrt(x);
        if ((r - Math.floor(r)) == 0) {
            return true;
        }
        return false;
    }
    
    
    static class LinkLNode{
        int data;
        LinkLNode next=null;
        LinkLNode(int x){
            this.data=x;
        }
    }
    static LinkLNode head = null;
        public static void insert_node(int position, int value) {
            // @params position, integer
            // @params value, integer
            int count = 0;
            LinkLNode trav = head;
            while(trav != null){
                trav=trav.next;
                count++;
            }
            LinkLNode newNode = new LinkLNode(value);
            if(position == 1){
                newNode.next = head;
                head = newNode;
            }else if(position <= count){
//            	int toprev = 1;
            	position--;
                trav = head;
                // going to 1 less node
                while( position > 1 ){
                    trav = trav.next;
                    position--;
//                    toprev++;
                }
            	LinkLNode nextNode = trav.next;
                // storing next node
                // adding next node
                trav.next = newNode;
                newNode.next = nextNode; 
            }else{
//            	last position
            	if(position == count+1) {
            		trav = head;
                    while(trav.next != null){
                        trav = trav.next;
                    }
                    trav.next = newNode;
            	}
                
            }
            
        }

        public static void delete_node(int position) {
            // @params position, integer
            int count=0;
            LinkLNode trav = head;
            while(trav != null){
                trav=trav.next;
                count++;
            }
            trav = head;
            //present positions of nodes            
            if(position == 1){
                head=trav.next;
            }else if(position < count){
                //as we already at head
            	position--;
                while(position > 1){
                    trav=trav.next;
                    position--;
                }
                // trav is 1 less than the actual node to del
                // storing noxt node, and to nullify toremove node, toremove.next=null
                LinkLNode nextNode = trav.next.next;
                trav.next = nextNode;

            }else{
            	//2nd last element
            	if(position == count) {
	                while(trav.next.next != null){
	                    trav=trav.next;
	                }
	                trav.next=null;
            	}
                
            }
            
        }

        public static void print_ll() {
            // Output each element followed by a space
            LinkLNode trav = head;
            while(trav != null){
                System.out.print(trav.data+" ");
                trav=trav.next;
            }
            System.out.println();
        }
    
public static void main(String[] args) {
        ArraySol<Integer> sol = new ArraySolClass<Integer>();
        insert_node(1, 1);
        insert_node(2, 2);
        insert_node(3, 3);
        insert_node(4, 4);
        print_ll();
        delete_node(4);
        print_ll();
        delete_node(1);
        print_ll();
        
        
        //int[] a = { 1, 1, 2, 2 };
    	//int[] b = { 1, 2, 1, 2 };
//        ArrayList<Integer> rest1 = new ArrayList<>(Arrays.asList( 16777, 1179, 265, 135, 90, 135, 34));
//        int[] rest = { 8, 8, 10, 5, 9, 1, 10, 4, 3, 5, 7 };
        //        
//        ArrayList<Integer> rest1 = new ArrayList<>(Arrays.asList( 4, 5, 2, 9, 8, 3, 10));
        
//        showAL(smallerIndexFromLeft);
////        showAL(smallerIndexFromRight);
////        showAL(greaterIndexFromLeft);
//        showAL(greaterIndexFromRight);
        
//        final ArrayList<Integer> rest1 = new ArrayList<>(Arrays.asList(  -49, 33, 35, 42   ));
//        
//        final ArrayList<Integer> rest2 = new ArrayList<>(Arrays.asList( -26  ));
//        String res = "(a+b-c-d+e-f+g+h+i)";
//        String giv = "xxikrwmjvsvckfrqxnibkcasompsuyuogauacjrr";
//        String giv1 = "wtcyawv";
//        ArrayList<ArrayList<Integer>> A = new ArrayList<ArrayList<Integer>> ();
//        A.add(new ArrayList<>(Arrays.asList( 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3)));
//        A.add(new ArrayList<>(Arrays.asList( 2, 6, 9)));
//        A.add(new ArrayList<>(Arrays.asList( 3, 6, 9)));
//        A.add(new ArrayList<>(Arrays.asList( 2, 0)));
//        A.add(new ArrayList<>(Arrays.asList( 1, 6)));
//        A.add(new ArrayList<>(Arrays.asList( 1, 6)));
//        A.add(new ArrayList<>(Arrays.asList( 2, 0)));
//        A.add(new ArrayList<>(Arrays.asList( 1, 5)));
//        A.add(new ArrayList<>(Arrays.asList( 1, 1)));
//        
//        StrVal[] newClassArr = new StrVal[4];
//        newClassArr[0] = new StrVal(9);
//        newClassArr[1] = new StrVal(123);
//        newClassArr[2] = new StrVal(45);
//        newClassArr[3] = new StrVal(2);
//        Arrays.sort(newClassArr);
//        for(int i=0;i < newClassArr.length;i++) {
//        	System.out.print(newClassArr[i].number+ " , ");
//        }
        
//        System.out.println( findMedian(A));
        
        
//        System.out.println(solve(giv));
        //it must be my issue, i did not get it the resolve.. again got same TA(was hoping new one 2nd time)
        
        
        
    }

    
}

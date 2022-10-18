/*
This program finds all occurences of a pattern in a text
It uses three algorithms and displays the time for each algorithm
By: Nividh Singh
Date: 10/18/2022

Most of the code for the searching was taken from geeksforgeeks and modified to work in this program
*/

// A prime number used by the Robin-Karp Algorithm
let q = 101;
// d is the number of characters in the input alphabet used by the the Z Algorithm
let d = 256;

// Code for if the button is clicked
document.getElementById("find").onclick = buttonClicked;

// Function for if the button is clicked
function buttonClicked() {

	// Gets the pattern and text from input boxes
  var pattern = document.getElementById("p").value;
  var text = document.getElementById("t").value;
  
  // If the pattern and text length is greater than 0 then we search
  if (pattern.length > 0 && text.length > 0) {
  	
    // variable for output field
    var outputField = document.getElementById("output")
    
    // Prints the results of search
    outputField.innerHTML = "The pattern is present at the following indexes: " + RKSearch(pattern, text, q);
    
    // Finds time for each algorithm and prints it
    const RKstart = Date.now();
    RKSearch(pattern, text, q);
    const RKend = Date.now();
    outputField.innerHTML += "<p></p>" + "The Robin Karp Algorithm took " + (RKend - RKstart) + " ms"
    const Zstart = Date.now();
    ZSearch(pattern, text);
    const Zend = Date.now();
    outputField.innerHTML += "<p></p>" + "The Z Algorithm took " + (Zend - Zstart) + " ms"
    const KMPstart = Date.now();
    KMPSearch(pattern, text);
    const KMPend = Date.now();
    outputField.innerHTML += "<p></p>" + "The KMP Algorithm took " + (KMPend - KMPstart) + " ms"

  }
}

// KMP algorithm based on https://www.geeksforgeeks.org/kmp-algorithm-for-pattern-searching/
//Finds the LPS Array
function computeLPSArray(pat, M, lps) {
  // length of the previous longest prefix suffix
  var len = 0;
  var i = 1;
  lps[0] = 0; // lps[0] is always 0

  // the loop calculates lps[i] for i = 1 to M-1
  while (i < M) {
    if (pat.charAt(i) == pat.charAt(len)) {
      len++;
      lps[i] = len;
      i++;
    } else // (pat[i] != pat[len])
    {
      // This is tricky. Consider the example.
      // AAACAAAA and i = 7. The idea is similar
      // to search step.
      if (len != 0) {
        len = lps[len - 1];

        // Also, note that we do not increment
        // i here
      } else // if (len == 0)
      {
        lps[i] = len;
        i++;
      }
    }
  }
}

//Searches using the LPS Array
function KMPSearch(pat, txt) {
  const indexes = []
  var M = pat.length;
  var N = txt.length;

  // create lps[] that will hold the longest
  // prefix suffix values for pattern
  var lps = [];
  var j = 0; // index for pat[]

  // Preprocess the pattern (calculate lps[]
  // array)
  computeLPSArray(pat, M, lps);

  var i = 0; // index for txt[]
  while ((N - i) >= (M - j)) {
    if (pat.charAt(j) == txt.charAt(i)) {
      j++;
      i++;
    }
    if (j == M) {
      indexes.push(i - j)
      //document.write("Found pattern " + "at index " + (i - j) + "\n");
      j = lps[j - 1];
    }

    // mismatch after j matches
    else if (i < N && pat.charAt(j) != txt.charAt(i)) {
      // Do not match lps[0..lps[j-1]] characters,
      // they will match anyway
      if (j != 0)
        j = lps[j - 1];
      else
        i = i + 1;
    }
  }
}


// Rabin Karp Algorithm based on https://www.geeksforgeeks.org/rabin-karp-algorithm-for-pattern-searching/ 
function RKSearch(pat, txt, q) {
  const indexes = [];
  let M = pat.length;
  let N = txt.length;
  let i, j;

  // Hash value for pattern
  let p = 0;

  // Hash value for txt
  let t = 0;
  let h = 1;

  // The value of h would be "pow(d, M-1)%q"
  for (i = 0; i < M - 1; i++)
    h = (h * d) % q;

  // Calculate the hash value of pattern and
  // first window of text
  for (i = 0; i < M; i++) {
    p = (d * p + pat[i].charCodeAt()) % q;
    t = (d * t + txt[i].charCodeAt()) % q;
  }

  // Slide the pattern over text one by one
  for (i = 0; i <= N - M; i++) {

    // Check the hash values of current
    // window of text and pattern. If the
    // hash values match then only
    // check for characters one by one
    if (p == t) {

      /* Check for characters one by one */
      for (j = 0; j < M; j++) {
        if (txt[i + j] != pat[j])
          break;
      }

      // if p == t and pat[0...M-1] =
      // txt[i, i+1, ...i+M-1]
      if (j == M)
        indexes.push(i)
      //document.write("Pattern found at index " + i + "<br/>");
    }

    // Calculate hash value for next window
    // of text: Remove leading digit, add
    // trailing digit
    if (i < N - M) {
      t = (d * (t - txt[i].charCodeAt() * h) +
        txt[i + M].charCodeAt()) % q;

      // We might get negative value of t,
      // converting it to positive
      if (t < 0)
        t = (t + q);
    }
  }
  return indexes;
}


//Z Algorithm based on https://www.geeksforgeeks.org/z-algorithm-linear-time-pattern-searching-algorithm/
//Searches using Z array
function ZSearch(pattern, text) {
  const indexes = []
  // Create concatenated string "P$T"
  let concat = pattern + "$" + text;

  let l = concat.length;

  let Z = new Array(l);

  // Construct Z array
  getZarr(concat, Z);

  // now looping through Z array for matching condition
  for (let i = 0; i < l; ++i) {

    // if Z[i] (matched region) is equal to pattern
    // length we got the pattern

    if (Z[i] == pattern.length) {
      indexes.push(i - pattern.length - 1)
      //document.write("Pattern found at index " + (i - pattern.length - 1) + "<br>");
      //alert("Pattern found at index " + (i - pattern.length - 1))
    }
  }
  return indexes
}

// Fills Z array for given string str[]
function getZarr(str, Z) {
  let n = str.length;

  // [L,R] make a window which matches with
  // prefix of s
  let L = 0,
    R = 0;

  for (let i = 1; i < n; ++i) {

    // if i>R nothing matches so we will calculate.
    // Z[i] using naive way.
    if (i > R) {

      L = R = i;

      // R-L = 0 in starting, so it will start
      // checking from 0'th index. For example,
      // for "ababab" and i = 1, the value of R
      // remains 0 and Z[i] becomes 0. For string
      // "aaaaaa" and i = 1, Z[i] and R become 5

      while (R < n && str[R - L] == str[R])
        R++;

      Z[i] = R - L;
      R--;

    } else {

      // k = i-L so k corresponds to number which
      // matches in [L,R] interval.
      let k = i - L;

      // if Z[k] is less than remaining interval
      // then Z[i] will be equal to Z[k].
      // For example, str = "ababab", i = 3, R = 5
      // and L = 2
      if (Z[k] < R - i + 1)
        Z[i] = Z[k];

      // For example str = "aaaaaa" and i = 2, R is 5,
      // L is 0
      else {


        // else start from R and check manually
        L = i;
        while (R < n && str[R - L] == str[R])
          R++;

        Z[i] = R - L;
        R--;
      }
    }
  }
}

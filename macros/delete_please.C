using namespace std;

void delete_please () {
  const int size = 5;
  
  int my_ints[size] = {1,3,2,5,7};
  set<int> test (my_ints, my_ints + size);

  set<int>::iterator it;
  for( it = test.begin(); it!=test.end(); ++it){
    int ans = *it;
    cout << ans << endl;
  }
  it = test.end();
  it --;
  int testit = *it; cout << testit << endl;
  return;
}

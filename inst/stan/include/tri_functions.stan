
  matrix to_lower_tri(vector x, int nr, int nc){
    matrix[nr,nc] y;
    int pos = 1;
    for(i in 1:nr){
      for(j in 1:nc){
        if(i < j){
          y[i,j] = 0;
        } else{
          y[i,j] = x[pos];
          pos += 1;
        }
      }
    }
    return y;
  }

  vector lt_to_vector(matrix y){
    int nr = rows(y);
    int nc = cols(y);
    vector[nc*(nr-nc) + min(rows(y),cols(y)) * (min(rows(y),cols(y)) + 1) / 2] x;
    int pos = 1;
    // now fill vector with matrix elements
    for(i in 1:nr){
      for(j in 1:nc){
        if(i < j) continue;

        x[pos] = y[i,j];
        pos += 1;

      }
    }
    return x;
  }


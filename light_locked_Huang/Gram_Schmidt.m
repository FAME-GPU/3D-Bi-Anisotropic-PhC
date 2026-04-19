function [n_Q, Y] = Gram_Schmidt(n_Y, Y)

   %
   %  Reorthogonal 
   %        modified G.M.
   %

   n_Q = 1;
   for k = 1:n_Y
      temp = norm( Y(:,k) );
      if ( abs(temp) > 1.0e-6 )
         Y(:,n_Q) = Y(:,k) / temp;
         if ( k < n_Y )
            for j = k+1:n_Y
               temp   = Y(:,n_Q)' * Y(:,j);
               Y(:,j) = Y(:,j) - temp * Y(:,n_Q);
            end
         end
         n_Q = n_Q + 1;
      end
   end
   n_Q = n_Q - 1;

   Y   = Y(:,1:n_Q);


% https://github.com/tal-amir/sparse-approximation-gsm
% sinteractive -p broadwl --mem=24G -c 4 --time=80:00:00
% module load gcc/10.2.0
% module load matlab/2022b
% matlab -nosplash -nodesktop
addpath('trimmed_lasso')
addpath('trimmed_lasso/gsm')
for ss = [1 5 20 100 500 2000 10000]
  rng(1);
  for ii = 1:20
    nn = "matdata/n500p10000s" + ss + "_normal_seed" + ii + ".mat";
    load(nn);
    disp("file:" + nn);
    disp("sparsity:" + ss);
    disp("ii:" + ii);
    k = [1 2 3 5 10 20 100 500];
    n = length(k);
    m = size(X,2);
    cv = zeros(n,1);
    B = zeros(m,n);
    t = 0;
    aa = randperm(500);
    for i = 1:n
      disp("k:" + k(i));
      cv_error = 0;
      X1 = X(aa,:); y1 = y(aa);
      X2 = X1(1:100,:); y2 = y1(1:100);
      X1(1:100,:) = [];
      y1(1:100) = [];
      evalc("[b1 out] = sparse_approx_gsm_v1_22(X1,y1,k(i),'profile','fast')");
      cv_error = cv_error + norm(X2*b1 - y2)^2;
      t = t + out.tElapsed;
      
      X1 = X(aa,:); y1 = y(aa);
      X2 = X1(101:200,:); y2 = y1(101:200);
      X1(101:200,:) = [];
      y1(101:200) = [];
      evalc("[b1 out] = sparse_approx_gsm_v1_22(X1,y1,k(i),'profile','fast')");
      cv_error = cv_error + norm(X2*b1 - y2)^2;
      t = t + out.tElapsed;

      X1 = X(aa,:); y1 = y(aa);
      X2 = X1(201:300,:); y2 = y1(201:300);
      X1(201:300,:) = [];
      y1(201:300) = [];
      evalc("[b1 out] = sparse_approx_gsm_v1_22(X1,y1,k(i),'profile','fast')");
      cv_error = cv_error + norm(X2*b1 - y2)^2;
      t = t + out.tElapsed;

      X1 = X(aa,:); y1 = y(aa);
      X2 = X1(301:400,:); y2 = y1(301:400);
      X1(301:400,:) = [];
      y1(301:400) = [];
      evalc("[b1 out] = sparse_approx_gsm_v1_22(X1,y1,k(i),'profile','fast')");
      cv_error = cv_error + norm(X2*b1 - y2)^2;
      t = t + out.tElapsed;

      X1 = X(aa,:); y1 = y(aa);
      X2 = X1(401:500,:); y2 = y1(401:500);
      X1(401:500,:) = [];
      y1(401:500) = [];
      evalc("[b1 out] = sparse_approx_gsm_v1_22(X1,y1,k(i),'profile','fast')");
      cv_error = cv_error + norm(X2*b1 - y2)^2;
      t = t + out.tElapsed;

      cv(i) = cv_error;
    end
    % disp(cv);
    [minval, minind] = min(cv);
    disp("minimum k:" + k(minind));
    evalc("[b out] = sparse_approx_gsm_v1_22(X,y,k(minind),'profile','fast')");
    t = t + out.tElapsed;
    disp("time:" + t);

    save("result_fast_n500p10000s" + ss + "_normal_seed" + ii + ".mat",'b');
    save("time_fast_n500p10000s" + ss + "_normal_seed" + ii + ".mat",'t');

    evalc("[b out] = sparse_approx_gsm_v1_22(X,y,20,'profile','fast')");
    tt = out.tElapsed;
    save("fixresult_fast_n500p10000s" + ss + "_normal_seed" + ii + ".mat",'b');
    save("fixtime_fast_n500p10000s" + ss + "_normal_seed" + ii + ".mat",'tt');
  end
end
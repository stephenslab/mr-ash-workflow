for p = [20 50 100 200 500 1000 2000 5000 10000 20000]
    rng(1);
    for ii = 1:20
        nn = "n500p" + p + "s20_largenfixedp_seed" + ii + ".mat";
        load(nn);
        disp("file:" + nn);
        disp("sparsity:" + ss);
        if p == 20
            k = [1 2 5 10 15 19];
        else
            k = [1 2 5 10 15 20 25 30];
        end
        n = length(k);
        m = size(X,2);
        cv = zeros(n,1);
        B = zeros(m,n);
        t = 0;
        aa = randperm(500);
        for i = 1:n
            disp(i);
            cv_error = 0;
            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(1:100,:); y2 = y1(1:100);
            X1(1:100,:) = [];
            y1(1:100) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(101:200,:); y2 = y1(101:200);
            X1(101:200,:) = [];
            y1(101:200) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(201:300,:); y2 = y1(201:300);
            X1(201:300,:) = [];
            y1(201:300) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(301:400,:); y2 = y1(301:400);
            X1(301:400,:) = [];
            y1(301:400) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(401:500,:); y2 = y1(401:500);
            X1(401:500,:) = [];
            y1(401:500) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            cv(i) = cv_error;
        end
        disp(cv);
        [minval, minind] = min(cv);
        disp("minimum k:" + k(minind));
        evalc("[b out] = sparse_approx_gsm(X,y,k(minind),'profile','fast')");
        t = t + out.tElapsed;
        disp("time:" + t);

        save("result_fast_n500p" + p + "s20_largenfixedp_seed" + ii + ".mat", 'b');
        save("time_fast_n500p" + p + "s20_largenfixedp_seed" + ii + ".mat", 't');
    end
end

for ss = [1 5 20 100 500]
    rng(1);
    for ii = 1:20
        nn = "n287s" + ss + "_geno_seed" + ii + ".mat";
        load(nn);
        disp("file:" + nn);
        disp("sparsity:" + ss);
        k = [1 2 5 10 20 30 50 100];
        n = length(k);
        m = size(X,2);
        cv = zeros(n,1);
        B = zeros(m,n);
        t = 0;
        aa = randperm(287);
        for i = 1:n
            disp(i);
            cv_error = 0;
            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(1:58,:); y2 = y1(1:58);
            X1(1:58,:) = [];
            y1(1:58) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(59:116,:); y2 = y1(59:116);
            X1(59:116,:) = [];
            y1(59:116) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(117:174,:); y2 = y1(117:174);
            X1(117:174,:) = [];
            y1(117:174) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(175:231,:); y2 = y1(175:231);
            X1(175:231,:) = [];
            y1(175:231) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(232:287,:); y2 = y1(232:287);
            X1(232:287,:) = [];
            y1(232:287) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            cv(i) = cv_error;
        end
        disp(cv);
        [minval, minind] = min(cv);
        disp("minimum k:" + k(minind));
        evalc("[b out] = sparse_approx_gsm(X,y,k(minind),'profile','fast')");
        t = t + out.tElapsed;
        disp("time:" + t);

        save("result_fast_n287s" + ss + "_geno_seed" + ii + ".mat", 'b');
        save("time_fast_n287s" + ss + "_geno_seed" + ii + ".mat", 't');

        evalc("[b out] = sparse_approx_gsm(X,y,20,'profile','fast')");
        tt = out.tElapsed;
        save("fixresult_fast_n287s" + ss + "_geno_seed" + ii + ".mat", 'b');
        save("fixtime_fast_n287s" + ss + "_geno_seed" + ii + ".mat", 'tt');
    end
end


for ss = [1 2 5 10 20 50 100 200]
    rng(1);
    for ii = 1:20
        nn = "n500p200s" + ss + "_normal_seed" + ii + ".mat";
        load(nn);
        disp("file:" + nn);
        disp("sparsity:" + ss);
        k = [1 2 5 10 20 30 50 100];
        n = length(k);
        m = size(X,2);
        cv = zeros(n,1);
        B = zeros(m,n);
        t = 0;
        aa = randperm(500);
        for i = 1:n
            disp(i);
            cv_error = 0;
            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(1:100,:); y2 = y1(1:100);
            X1(1:100,:) = [];
            y1(1:100) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(101:200,:); y2 = y1(101:200);
            X1(101:200,:) = [];
            y1(101:200) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(201:300,:); y2 = y1(201:300);
            X1(201:300,:) = [];
            y1(201:300) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(301:400,:); y2 = y1(301:400);
            X1(301:400,:) = [];
            y1(301:400) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(401:500,:); y2 = y1(401:500);
            X1(401:500,:) = [];
            y1(401:500) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            cv(i) = cv_error;
        end
        disp(cv);
        [minval, minind] = min(cv);
        disp("minimum k:" + k(minind));
        evalc("[b out] = sparse_approx_gsm(X,y,k(minind),'profile','fast')");
        t = t + out.tElapsed;
        disp("time:" + t);

        save("result_fast_n500p200s" + ss + "_normal_seed" + ii + ".mat", 'b');
        save("time_fast_n500p200s" + ss + "_normal_seed" + ii + ".mat", 't');

        evalc("[b out] = sparse_approx_gsm(X,y,20,'profile','fast')");
        tt = out.tElapsed;
        save("fixresult_fast_n500p200s" + ss + "_normal_seed" + ii + ".mat", 'b');
        save("fixtime_fast_n500p200s" + ss + "_normal_seed" + ii + ".mat", 'tt');
    end
end
 
for ss = [1 3 10 30 100 300 1000]
    rng(1);
    for ii = 1:20
        nn = "n500p1000s" + ss + "_normal_seed" + ii + ".mat";
        load(nn);
        disp("file:" + nn);
        disp("sparsity:" + ss);
        k = [1 2 3 5 10 20 30];
        n = length(k);
        m = size(X,2);
        cv = zeros(n,1);
        B = zeros(m,n);
        t = 0;
        aa = randperm(500);
        for i = 1:n
            disp(i);
            cv_error = 0;
            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(1:100,:); y2 = y1(1:100);
            X1(1:100,:) = [];
            y1(1:100) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(101:200,:); y2 = y1(101:200);
            X1(101:200,:) = [];
            y1(101:200) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(201:300,:); y2 = y1(201:300);
            X1(201:300,:) = [];
            y1(201:300) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(301:400,:); y2 = y1(301:400);
            X1(301:400,:) = [];
            y1(301:400) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(401:500,:); y2 = y1(401:500);
            X1(401:500,:) = [];
            y1(401:500) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            cv(i) = cv_error;
        end
        disp(cv);
        [minval, minind] = min(cv);
        disp("minimum k:" + k(minind));
        evalc("[b out] = sparse_approx_gsm(X,y,k(minind),'profile','fast')");
        t = t + out.tElapsed;
        disp("time:" + t);

        save("result_fast_n500p1000s" + ss + "_normal_seed" + ii + ".mat", 'b');
        save("time_fast_n500p1000s" + ss + "_normal_seed" + ii + ".mat", 't');

        evalc("[b out] = sparse_approx_gsm(X,y,20,'profile','fast')");
        tt = out.tElapsed;
        save("fixresult_fast_n500p1000s" + ss + "_normal_seed" + ii + ".mat", 'b');
        save("fixtime_fast_n500p1000s" + ss + "_normal_seed" + ii + ".mat", 'tt');
    end
end

for ss = [1 5 20 100 500 2000 10000]
    rng(1);
    for ii = 1:20
        nn = "n500p10000s" + ss + "_normal_seed" + ii + ".mat";
        load(nn);
        disp("file:" + nn);
        disp("sparsity:" + ss);
        k = [1 2 3 5 10 20 100 500];
        n = length(k);
        m = size(X,2);
        cv = zeros(n,1);
        B = zeros(m,n);
        t = 0;
        aa = randperm(500);
        for i = 1:n
            disp(i);
            cv_error = 0;
            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(1:100,:); y2 = y1(1:100);
            X1(1:100,:) = [];
            y1(1:100) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(101:200,:); y2 = y1(101:200);
            X1(101:200,:) = [];
            y1(101:200) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(201:300,:); y2 = y1(201:300);
            X1(201:300,:) = [];
            y1(201:300) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(301:400,:); y2 = y1(301:400);
            X1(301:400,:) = [];
            y1(301:400) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(401:500,:); y2 = y1(401:500);
            X1(401:500,:) = [];
            y1(401:500) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            cv(i) = cv_error;
        end
        disp(cv);
        [minval, minind] = min(cv);
        disp("minimum k:" + k(minind));
        evalc("[b out] = sparse_approx_gsm(X,y,k(minind),'profile','fast')");
        t = t + out.tElapsed;
        disp("time:" + t);

        save("result_fast_n500p10000s" + ss + "_normal_seed" + ii + ".mat", 'b');
        save("time_fast_n500p10000s" + ss + "_normal_seed" + ii + ".mat", 't');

        evalc("[b out] = sparse_approx_gsm(X,y,20,'profile','fast')");
        tt = out.tElapsed;
        save("fixresult_fast_n500p10000s" + ss + "_normal_seed" + ii + ".mat", 'b');
        save("fixtime_fast_n500p10000s" + ss + "_normal_seed" + ii + ".mat", 'tt');
    end
end


for ss = [1 3 10 30 100 300 1000]
    rng(1);
    for ii = 1:20
        nn = "n500p1000s" + ss + "_rho095_seed" + ii + ".mat";
        load(nn);
        disp("file:" + nn);
        disp("sparsity:" + ss);
        k = [1 2 3 5 10 20 30];
        n = length(k);
        m = size(X,2);
        cv = zeros(n,1);
        B = zeros(m,n);
        t = 0;
        aa = randperm(500);
        for i = 1:n
            disp(i);
            cv_error = 0;
            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(1:100,:); y2 = y1(1:100);
            X1(1:100,:) = [];
            y1(1:100) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(101:200,:); y2 = y1(101:200);
            X1(101:200,:) = [];
            y1(101:200) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(201:300,:); y2 = y1(201:300);
            X1(201:300,:) = [];
            y1(201:300) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(301:400,:); y2 = y1(301:400);
            X1(301:400,:) = [];
            y1(301:400) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(401:500,:); y2 = y1(401:500);
            X1(401:500,:) = [];
            y1(401:500) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            cv(i) = cv_error;
        end
        disp(cv);
        [minval, minind] = min(cv);
        disp("minimum k:" + k(minind));
        evalc("[b out] = sparse_approx_gsm(X,y,k(minind),'profile','fast')");
        t = t + out.tElapsed;
        disp("time:" + t);

        save("result_fast_n500p1000s" + ss + "_rho095_seed" + ii + ".mat", 'b');
        save("time_fast_n500p1000s" + ss + "_rho095_seed" + ii + ".mat", 't');

        evalc("[b out] = sparse_approx_gsm(X,y,20,'profile','fast')");
        tt = out.tElapsed;
        save("fixresult_fast_n500p1000s" + ss + "_rho095_seed" + ii + ".mat", 'b');
        save("fixtime_fast_n500p1000s" + ss + "_rho095_seed" + ii + ".mat", 'tt');
    end
end

for ss = [1 3 10 30 100 300 1000]
    rng(1);
    for ii = 1:20
        nn = "n500p1000s" + ss + "_rho050_seed" + ii + ".mat";
        load(nn);
        disp("file:" + nn);
        disp("sparsity:" + ss);
        k = [1 2 3 5 10 20 30];
        n = length(k);
        m = size(X,2);
        cv = zeros(n,1);
        B = zeros(m,n);
        t = 0;
        aa = randperm(500);
        for i = 1:n
            disp(i);
            cv_error = 0;
            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(1:100,:); y2 = y1(1:100);
            X1(1:100,:) = [];
            y1(1:100) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(101:200,:); y2 = y1(101:200);
            X1(101:200,:) = [];
            y1(101:200) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(201:300,:); y2 = y1(201:300);
            X1(201:300,:) = [];
            y1(201:300) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(301:400,:); y2 = y1(301:400);
            X1(301:400,:) = [];
            y1(301:400) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            X1 = X(aa,:); y1 = y(aa);
            X2 = X1(401:500,:); y2 = y1(401:500);
            X1(401:500,:) = [];
            y1(401:500) = [];
            evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
            cv_error = cv_error + norm(X2*b1 - y2)^2;
            t = t + out.tElapsed;

            cv(i) = cv_error;
        end
        disp(cv);
        [minval, minind] = min(cv);
        disp("minimum k:" + k(minind));
        evalc("[b out] = sparse_approx_gsm(X,y,k(minind),'profile','fast')");
        t = t + out.tElapsed;
        disp("time:" + t);

        save("result_fast_n500p1000s" + ss + "_rho050_seed" + ii + ".mat", 'b');
        save("time_fast_n500p1000s" + ss + "_rho050_seed" + ii + ".mat", 't');

        evalc("[b out] = sparse_approx_gsm(X,y,20,'profile','fast')");
        tt = out.tElapsed;
        save("fixresult_fast_n500p1000s" + ss + "_rho050_seed" + ii + ".mat", 'b');
        save("fixtime_fast_n500p1000s" + ss + "_rho050_seed" + ii + ".mat", 'tt');
    end
end

for sett = ["t1", "t2", "t3", "t4", "lap", "normal", "unif", "const", "pve0.000000", "pve0.010000", "pve0.050000", "pve0.100000", "pve0.200000", "pve0.300000", "pve0.400000", "pve0.500000", "pve0.600000", "pve0.700000", "pve0.800000", "pve0.900000", "pve0.950000", "pve0.990000"]
    for ss = [1000]
        rng(1);
        for ii = 1:20
            nn = "n500p1000s" + ss + "_" + sett + "_seed" + ii + ".mat";
            load(nn);
            disp("file:" + nn);
            disp("sparsity:" + ss);
            k = [1 2 3 5 10 20 30];
            n = length(k);
            m = size(X,2);
            cv = zeros(n,1);
            B = zeros(m,n);
            t = 0;
            aa = randperm(500);
            for i = 1:n
                disp(i);
                cv_error = 0;
                X1 = X(aa,:); y1 = y(aa);
                X2 = X1(1:100,:); y2 = y1(1:100);
                X1(1:100,:) = [];
                y1(1:100) = [];
                evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
                cv_error = cv_error + norm(X2*b1 - y2)^2;
                t = t + out.tElapsed;
    
                X1 = X(aa,:); y1 = y(aa);
                X2 = X1(101:200,:); y2 = y1(101:200);
                X1(101:200,:) = [];
                y1(101:200) = [];
                evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
                cv_error = cv_error + norm(X2*b1 - y2)^2;
                t = t + out.tElapsed;
    
                X1 = X(aa,:); y1 = y(aa);
                X2 = X1(201:300,:); y2 = y1(201:300);
                X1(201:300,:) = [];
                y1(201:300) = [];
                evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
                cv_error = cv_error + norm(X2*b1 - y2)^2;
                t = t + out.tElapsed;
    
                X1 = X(aa,:); y1 = y(aa);
                X2 = X1(301:400,:); y2 = y1(301:400);
                X1(301:400,:) = [];
                y1(301:400) = [];
                evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
                cv_error = cv_error + norm(X2*b1 - y2)^2;
                t = t + out.tElapsed;
    
                X1 = X(aa,:); y1 = y(aa);
                X2 = X1(401:500,:); y2 = y1(401:500);
                X1(401:500,:) = [];
                y1(401:500) = [];
                evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
                cv_error = cv_error + norm(X2*b1 - y2)^2;
                t = t + out.tElapsed;
    
                cv(i) = cv_error;
            end
            disp(cv);
            [minval, minind] = min(cv);
            disp("minimum k:" + k(minind));
            evalc("[b out] = sparse_approx_gsm(X,y,k(minind),'profile','fast')");
            t = t + out.tElapsed;
            disp("time:" + t);
    
            save("result_fast_n500p1000s" + ss + "_" + sett + "_seed" + ii + ".mat", 'b');
            save("time_fast_n500p1000s" + ss + "_" + sett + "_seed" + ii + ".mat", 't');
    
            evalc("[b out] = sparse_approx_gsm(X,y,20,'profile','fast')");
            tt = out.tElapsed;
            save("fixresult_fast_n500p1000s" + ss + "_" + sett + "_seed" + ii + ".mat", 'b');
            save("fixtime_fast_n500p1000s" + ss + "_" + sett + "_seed" + ii + ".mat", 'tt');
        end
    end
end

for sett = ["t1", "t2", "t4", "t8", "lap", "normal", "unif"]
    for ss = [20, 500]
        rng(1);
        for ii = 1:20
            nn = "n500p1000s" + ss + "_noise=" + sett + "_seed" + ii + ".mat";
            load(nn);
            disp("file:" + nn);
            disp("sparsity:" + ss);
            k = [1 2 3 5 10 19];
            n = length(k);
            m = size(X,2);
            cv = zeros(n,1); 
            B = zeros(m,n);
            t = 0;
            aa = randperm(500);
            for i = 1:n
                disp(i);
                cv_error = 0;
                X1 = X(aa,:); y1 = y(aa);
                X2 = X1(1:100,:); y2 = y1(1:100);
                X1(1:100,:) = [];
                y1(1:100) = [];
                evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
                cv_error = cv_error + norm(X2*b1 - y2)^2;
                t = t + out.tElapsed;
    
                X1 = X(aa,:); y1 = y(aa);
                X2 = X1(101:200,:); y2 = y1(101:200);
                X1(101:200,:) = [];
                y1(101:200) = [];
                evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
                cv_error = cv_error + norm(X2*b1 - y2)^2;
                t = t + out.tElapsed;
    
                X1 = X(aa,:); y1 = y(aa);
                X2 = X1(201:300,:); y2 = y1(201:300);
                X1(201:300,:) = [];
                y1(201:300) = [];
                evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
                cv_error = cv_error + norm(X2*b1 - y2)^2;
                t = t + out.tElapsed;
    
                X1 = X(aa,:); y1 = y(aa);
                X2 = X1(301:400,:); y2 = y1(301:400);
                X1(301:400,:) = [];
                y1(301:400) = [];
                evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
                cv_error = cv_error + norm(X2*b1 - y2)^2;
                t = t + out.tElapsed;
    
                X1 = X(aa,:); y1 = y(aa);
                X2 = X1(401:500,:); y2 = y1(401:500);
                X1(401:500,:) = [];
                y1(401:500) = [];
                evalc("[b1 out] = sparse_approx_gsm(X1,y1,k(i),'profile','fast')");
                cv_error = cv_error + norm(X2*b1 - y2)^2;
                t = t + out.tElapsed;
    
                cv(i) = cv_error;
            end
            disp(cv);
            [minval, minind] = min(cv);
            disp("minimum k:" + k(minind));
            evalc("[b out] = sparse_approx_gsm(X,y,k(minind),'profile','fast')");
            t = t + out.tElapsed;
            disp("time:" + t);
    
            save("result_fast_n500p1000s" + ss + "_noise" + sett + "_seed" + ii + ".mat", 'b');
            save("time_fast_n500p1000s" + ss + "_noise" + sett + "_seed" + ii + ".mat", 't');
    
            evalc("[b out] = sparse_approx_gsm(X,y,20,'profile','fast')");
            tt = out.tElapsed;
            save("fixresult_fast_n500p1000s" + ss + "_noise" + sett + "_seed" + ii + ".mat", 'b');
            save("fixtime_fast_n500p1000s" + ss + "_noise" + sett + "_seed" + ii + ".mat", 'tt');
        end
    end
end

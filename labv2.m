

eps = 0.001;



[X, map] = imread('x1.bmp');
Y = double(imread('y8.bmp'));

X = double(X);

A_ = Y * pinv(X);
Y_ = A_ * X;

A__ = Y * penmo(X);
Y__ = A__ * X;

A___ = Y * grevi(X);
Y___ = A__ * X;


Y___ = uint8(Y___);


imwrite(Y___, map, 'ogy1.bmp', 'bmp');

function M = penmo(X)
	[m , n] = size(X);
	dlt = 1.0;
	M = X.';
	cond = true;
	while cond
		it = (X.' * X + dlt * eye(n)) \ X.';
		cond = norm(it - M) <= eps;
		M = it;
		dlt = dlt / 2;
	end
end

function Y = grevi(X)
	[m, n] = size(X);
	E = eye(n);
	a1 = X(1, :).';
	Y = a1 / (a1.' * a1);
	for i = 2:m
		a = X(i, :).';
		Z = E - Y * X(1:(i - 1), :);
		R = Y * Y.';
		if a.' * Z * a > 0
			Y = Y - (Z * a * a.' * Y) / (a.' * Z * a);
			Y(:, i) = (Z * a) / (a.' * Z * a);
		else
			Y = Y - (R * a * a.' * Y) / (1 + a.' * R * a);
			Y(:, i) = (R * a) / (1 + a.' * R * a);
		end
	end
end

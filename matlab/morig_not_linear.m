%%
close all

%%
rng(2019)

%%
dataset = [repmat([2, 2], 200, 1) + randn(200, 2) * chol([1 .5; .5 2]);
           repmat([-2, 2], 200, 1) + randn(200, 2) * chol([1 .8; .8 1.6]);
           repmat([2, -2], 200, 1) + randn(200, 2) * chol([2 .4; .4 1]);
           repmat([-2, -2], 200, 1) + randn(200, 2) * chol([1 .5; .5 1]);
          ];

yes = [zeros(400, 1) + 1; zeros(400, 1) + 2];
no = [zeros(200, 1) + 1; zeros(200, 1) + 2; zeros(200, 1) + 1; zeros(200, 1) + 2];
colors = ["r", "b"];

dataset = dataset - mean(dataset, 1);
dataset = dataset ./ std(dataset, 1);

figure
sgtitle('Preview')
subplot(1, 2, 1)
scatter(dataset(yes == 1, 1), dataset(yes == 1, 2), 2)
hold on
scatter(dataset(yes == 2, 1), dataset(yes == 2, 2), 2)
pbaspect([1 1 1])
title('colored by yes - want better separation')

subplot(1, 2, 2)
scatter(dataset(no == 1, 1), dataset(no == 1, 2), 2)
hold on
scatter(dataset(no == 2, 1), dataset(no == 2, 2), 2)
pbaspect([1 1 1])
title('colored by no - want less separation')

%%
S = 0;
for i = 1:2
    temp = dataset(no ~= i, :) - mean(dataset(no == 1, :), 1);
    S = S + temp' * temp;
end
S = S / size(dataset, 1);

%%
D = mpdist(dataset, S);

%u = umap();
%u.metric = 'precomputed';
%R = u.fit(D);
R = tsne(D, [], 2);

%%

figure
sgtitle('Corrected')
subplot(1, 2, 1)
scatter(R(yes == 1, 1), R(yes == 1, 2), 2)
hold on
scatter(R(yes == 2, 1), R(yes == 2, 2), 2)
pbaspect([1 1 1])
title('colored by yes - want better separation')

subplot(1, 2, 2)
scatter(R(no == 1, 1), R(no == 1, 2), 2)
hold on
scatter(R(no == 2, 1), R(no == 2, 2), 2)
pbaspect([1 1 1])
title('colored by no - want less separation')

%%
D0 = mpdist(dataset, eye(size(dataset, 2)));

%u = umap();
%u.metric = 'precomputed';
%R = u.fit(D);
R0 = tsne(D0, [], 2);

figure
sgtitle('Original')
subplot(1, 2, 1)
scatter(R0(yes == 1, 1), R0(yes == 1, 2), 2)
hold on
scatter(R0(yes == 2, 1), R0(yes == 2, 2), 2)
pbaspect([1 1 1])
title('colored by yes - want better separation')

subplot(1, 2, 2)
scatter(R0(no == 1, 1), R0(no == 1, 2), 2)
hold on
scatter(R0(no == 2, 1), R0(no == 2, 2), 2)
pbaspect([1 1 1])
title('colored by no - want less separation')
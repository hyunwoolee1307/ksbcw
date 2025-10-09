%% mat 파일 만들기
clc; clear
%% data 불러오기
KODC = [];
for i=1:40
    KODC = vertcat(KODC,readtable(strcat('KODC_',string(1979+i),'.xls'),...
        'VariableNamingRule','preserve'));
end
clearvars i
%% 변수 저장
cruiseline = KODC.("정선");
station = KODC.("정점");
depth = KODC.("관측수심(m)");
years = year(datetime(KODC.("관측일시(KST)")));
months = month(datetime(KODC.("관측일시(KST)")));
days = day(datetime(KODC.("관측일시(KST)")));
temps = KODC.("수온(℃)");
salt = KODC.("염분(psu)");
oxygen = KODC.("용존산소(ml/L)");
temps_qc = KODC.(10);
salt_qc = KODC.(12);
oxygen_qc = KODC.(14);
%% 변수 타입 변환
%Char to String
cruiseline = convertCharsToStrings(cruiseline);
station = convertCharsToStrings(station);
depth = convertCharsToStrings(depth);
temps = convertCharsToStrings(temps);
salt = convertCharsToStrings(salt);
oxygen = convertCharsToStrings(oxygen);
temps_qc = convertCharsToStrings(temps_qc);
salt_qc = convertCharsToStrings(salt_qc);
oxygen_qc = convertCharsToStrings(oxygen_qc);
%% 수온, 염분, DO 종합하기
T = [cruiseline,station,depth,years,months,days,...
    temps,temps_qc,latitude,longitude];
T = double(T);
save("Temperature","T")
clearvars T
S = [cruiseline,station,depth,years,months,days,...
    salt,salt_qc];
S =  double(S);
save("Salt","S");
clearvars S
DO = [cruiseline,station,depth,years,months,days,...
    oxygen,oxygen_qc];
DO = double(DO);
save("Oxygen","DO");
clearvars DO
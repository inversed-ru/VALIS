{ 
Copyright (c) Peter Karpov 2010 - 2017.

Usage of the works is permitted provided that this instrument is retained with 
the works, so that any entity that uses the works is notified of this instrument.

DISCLAIMER: THE WORKS ARE WITHOUT WARRANTY.
}
{$IFDEF FPC} {$MODE DELPHI} {$ENDIF} {$APPTYPE CONSOLE}
program VALIS; //////////////////////////////////////////////////////////////////////
{
>> Version: 1.4

>> Description
   Implementation of VALIS (Vote-Allocating Immune System) classification algorithm.
   Refer to the original papers [1, 2] and the project homepage [3] for details.

>> Author
   Peter Karpov
   Email    : PeterKarpov@inversed.ru
   Homepage : inversed.ru
   GitHub   : inversed-ru
   Twitter  : @inversed_ru
   
>> References
   1. VALIS, a novel immune-inspired classification algorithm.
      Karpov P.
   2. VALIS: an immune-inspired classification algorithm.
      Karpov P., Tonda, A., Squillero, G.
   3. VALIS homepage
      http://inversed.ru/AIS.htm
   
>> ToDo
   - Implement k nearest binding
   - Read the settings from a config file
   ! Monitor more population stats: mean and SD of fitness, radius, NBound, ???
   ! Try fitness-proportional instead of rank-proportional selection
   ! Test on synthetic data with validation set
   ? Return best population based on accuracy and coverage
   ? Which coverage formulation to use?
   ? Use tie-breaks
   - Various voting modes
   ! Incremental updates
   - Alternative fitness measures (linear and log gains)
   - Should the binding functions be defined by antibody parameters?
   - Selective pressure parameter
   ! Data normalization
   ? Special infinite radius antibody to avoid zero coverage problems
   ? Use priors instead of pseudocounts
   - Try other antibodies: linear and quadratic shapes, decision trees
   - Problems to try:
     - Character recognition
     - Text classification, sentiment analysis
     - Time series classification
     - Text generation
   
>> Changelog
   1.4  : 2017.12.16 ~ InvLibs compatibility
                     + CalcWeightMatrix procedure
                     ~ Cleanup
   1.3  : 2017.02.09 + Crossvalidation
                     ~ Made statistical multiplier parameter continous
                     ~ Moved NClasses and NVars variables into TAntigens
                     ~ Replaced PopSize setting with Abundance ratio
                     + Coverage rank multiplier selection (repSumWRankMul)
                     + Alternative fitness calculation mode
                     + Alternative statistical multiplier
                     ~ In case of 0 bound antibodies, the one with min D / R is used
                       for classification
                     + Additive smoothing for class distributions
                     + Coverage calculation
   1.2  : 2017.01.03 * Bug in dual replacement
                     + More binding functions
   1.1  : 2016.12.16 + New replacement method
                     * Incorrect sharing factor calculation
   1.0  : 2016.12.15 ~ Rewritten from scratch
   0.0  : 2009.05.21 + Initial version
   Notation: + added, - removed, * fixed, ~ changed
}////////////////////////////////////////////////////////////////////////////////////

uses
      SysUtils,
      InvSys,
      Math,
      ExtraMath,
      Arrays,
      Sorting,
      StringUtils,
      RandVars,
      Statistics,
      Problem;

type
      TBindingFunction  =
         (bfSquare, bfTriangle, bfCircle, bfPoly, bfExp1, bfExp2, bfExp4);

      TTrainParams =
         record
         UseCrossover      :  Boolean;
         BindingFunction   :  TBindingFunction;
         MaxGens           :  Integer;         
         end;

      TReplacement = (repFitness, repDual, repSumWRankMul, repProdRank, repSumRank);

      TWeightMatrix = array of TRealArray;

const
      // True  2017 version, individual bond sharing factors
      // False 2009 version, average sharing factor
      AltFitness        =  True;
      
      BindingFunction   =  bfSquare;
      UseCrossover      =  True;
      RankProp          =  True;
      StatMul           =  0;
      AltStatMul        =  False;
      Pseudocount       =  1;
      Replacement       =  repFitness;
      SumWFrac          =  8;
      Abundance         =  1;
      MaxGens           =  600;
      SaveVis           =  False;
      PathVis           =  'Visualization' + PathDelim;
      DirData           =  'Datasets';

{-----------------------<< Binding >>-----------------------------------------------}

// Return the binding weight given the antibody-antigen distance D, 
// the antibody radius R, and the binding function F
function BindingWeight(
         D, R     :  Real;
         F        :  TBindingFunction
         )        :  Real;
   var
         x        :  Real;
   begin
   Assert(R > 0);
   x := D / R;
   case F of
      bfSquare    :  Result := UnitStep(1 - x);
      bfTriangle  :  Result := Max(1 - x, 0);
      bfCircle    :  Result := SafeSqrt(1 - Sqr(x));
      bfPoly      :  Result := (1 - 3 * x * x + 2 * x * x) * UnitStep(1 - x);
      bfExp1      :  Result := Exp(-x);
      bfExp2      :  Result := Exp(-Sqr(x));
      bfExp4      :  Result := Exp(-Sqr(Sqr(x)));
      else           Result := 0;
                     Assert(False);
      end;
   end;
   
   
// Statistical weight multiplier
function WeightMul(
      SumW  :  Real
      )     :  Real;
   begin
   if AltStatMul then
      Result := 1 - Exp(-SumW)
   else if StatMul = 0 then
      Result := 1
   else
      Result := Power(1 - Sqrt(1 / Max(SumW, 1)), StatMul);
   end;
   
   
// Calculate the antibody-antigen binding weights matrix W
procedure CalcWeightMatrix(
   var   W           :  TWeightMatrix;
   const Antibodies  :  TAntibodies;
   const Antigens    :  TAntigens);
   var
         i, j        :  Integer;
   begin
   SetLength(W, Antibodies.N, Antigens.N);
   for i := 0 to Antibodies.N - 1 do
      for j := 0 to Antigens.N - 1 do
         W[i, j] := BindingWeight(
            Distance(Antigens._[j], Antibodies._[i]),
            Antibodies._[i].Radius, BindingFunction);
   end;

{-----------------------<< Accuracy and antibody statistics >>----------------------}

// Return the classification accuracy given a weight matrix W
function CalcAccuracy(
   const Antibodies  :  TAntibodies;
   const Antigens    :  TAntigens;
   const W           :  TWeightMatrix
         )           :  Real;
   var
         i, j, 
         NCorrect,
         Voted       :  Integer;
         Votes       :  TRealArray;
         D, RelD,
         MinRelD     :  Real;
   begin
   NCorrect := 0;
   for j := 0 to Antigens.N - 1 do
      begin
      // Count the votes
      InitArray(Votes, Antigens.NClasses, 0);
      for i := 0 to Antibodies.N - 1 do
         if Antibodies._[i].SumW <> 0 then
            ArrayMulAdd(Votes, Antibodies._[i].ClassDistr, W[i, j]);

      // In case of no bound antibodies, use one with the min relative distance
      if Sum(Votes) = 0 then
         begin
         MinRelD := Infinity;
         for i := 0 to Antibodies.N - 1 do
            begin
            D := Distance(Antigens._[j], Antibodies._[i]);
            RelD := D / Antibodies._[i].Radius;
            if RelD < MinRelD then
               begin
               MinRelD := RelD;
               Votes := Copy(Antibodies._[i].ClassDistr);
               end;
            end;
         end;

      // Select the class with them most votes as the classification result
      Voted := RandMaxIndex(Votes);
      NCorrect := NCorrect + Ord(Antigens._[j].Class_ = Voted);
      end;
   Result := NCorrect / Antigens.N;
   end;


// Return the classification accuracy
function TestAccuracy(
   const Antibodies  :  TAntibodies;
   const Antigens    :  TAntigens
         )           :  Real;
   var
         W           :  TWeightMatrix;
   begin
   CalcWeightMatrix(W, Antibodies, Antigens);
   Result := CalcAccuracy(Antibodies, Antigens, W); 
   end;

(*
// Update Antibody - Antigens interaction statistics
// #TODO implement
procedure UpdateStats(
   var   Antibody :  TAntibody;
   const Antigens :  TAntigens);
   begin

   end;
*)

// Calculate the antibody statistics
// #HACK Huge
procedure CalcStats(
   var   Antibodies     :  TAntibodies;
   const Antigens       :  TAntigens);
   var
         i, j, k        :  Integer;
         W              :  TWeightMatrix;
         SumWAG, F      :  TRealArray;
         RelW,
         BondAccuracy   :  Real;
   begin
   // Calculate total antigen binding weights and coverage
   CalcWeightMatrix(W, Antibodies, Antigens);
   with Antibodies do
      begin
      Coverage := 0;
      InitArray(SumWAG, Antigens.N, 0);
      for j := 0 to Antigens.N - 1 do
         begin
         for i := 0 to Antibodies.N - 1 do
            SumWAG[j] := SumWAG[j] + W[i, j];
         Coverage := Coverage + 1 - 1 / (1 + SumWAG[j]);
         end;
      Coverage := Coverage / Antigens.N;
      end;

   // Calculate antibody total weights and class distributions
   for i := 0 to Antibodies.N - 1 do
      with Antibodies._[i] do
         begin
         SumW := 0;
         InitArray(ClassDistr, Antigens.NClasses, Pseudocount);
         for j := 0 to Antigens.N - 1 do
            begin
            SumW := SumW + W[i, j];
            k := Antigens._[j].Class_;
            ClassDistr[k] := ClassDistr[k] + W[i, j];
            end;
         Normalize(ClassDistr, normUnitSum);
         MajorityClass := RandMaxIndex(ClassDistr);
         end;

   // Calculate antibody accuracy
   for i := 0 to Antibodies.N - 1 do
      with Antibodies._[i] do
         begin
         Accuracy := 0;
         for j := 0 to Antigens.N - 1 do
            Accuracy := Accuracy + ClassDistr[Antigens._[j].Class_] * W[i, j];
         Accuracy := SafeDiv(Accuracy, SumW, 0);
         end;

   // Calculate the sharing factors
   for i := 0 to Antibodies.N - 1 do
      with Antibodies._[i] do
         begin
         KShare := 0;
         for j := 0 to Antigens.N - 1 do
            begin
            Assert(SumWAG[j] >= W[i, j]);
            KShare := KShare + SafeDiv(Sqr(W[i, j]), SumWAG[j], 0);
            end;
         KShare := SafeDiv(KShare, SumW, 1);
         end;

   // Sort the antibodies by total binding weight
   SetLength(F, Antibodies.N);
   for i := 0 to Antibodies.N - 1 do
      F[i] := Antibodies._[i].SumW;
   OrderRealArray(Antibodies.OrderSumW, F, soAScending);

   // Calculate fitness
   for i := 0 to Antibodies.N - 1 do
      with Antibodies._[i] do
         if AltFitness then
            begin
            Fitness := 0;
            for j := 0 to Antigens.N - 1 do
               begin
               RelW := SafeDiv(W[i, j], SumWAG[j], 0);
               BondAccuracy := ClassDistr[Antigens._[j].Class_];
               Fitness := Fitness + BondAccuracy * RelW * W[i, j];
               end;
            Fitness := SafeDiv(Fitness * WeightMul(SumW), SumW, 0);
            end
         else
            Fitness := KShare * Accuracy * WeightMul(SumW);

   // Apply SumW rank multipliers if necessary
   if Replacement = repSumWRankMul then
      for i := 0 to Antibodies.N - 1 do
         with Antibodies._[Antibodies.OrderSumW[i]] do
            Fitness := Fitness * (1 - Power(1 - i / (Antibodies.N - 1), SumWFrac));

   // Sort the antibodies by fitness
   SetLength(F, Antibodies.N);
   for i := 0 to Antibodies.N - 1 do
      F[i] := Antibodies._[i].Fitness;
   OrderRealArray(Antibodies.OrderFitness, F, soAScending);

   // Set the overall classification accuracy
   Antibodies.Accuracy := CalcAccuracy(Antibodies, Antigens, W);
   end;

{-----------------------<< Training >>----------------------------------------------}

// Return an index of a random antibody,
// use rank proportional selection if RankProp is on
function RandAntibodyIndex(
   const Antibodies  :  TAntibodies;
         RankProp    :  Boolean
         )           :  Integer;
   var
         i           :  Integer;
   begin
   with Antibodies do
      begin
      if RankProp then
         i := Max(Random(N), Random(N)) else
         i :=     Random(N);
      Result := OrderFitness[i];
      end;
   end;


// Create Children from Parents
procedure Replicate(
   var   Children       :  TAntibodies;
   const Parents        :  TAntibodies;
         NChildren      :  Integer;
         MutationRate   :  Real;
         UseCrossover   :  Boolean);
   var
         IdParentA,
         IdParentB, i   :  Integer;
   begin
   Children.N := NChildren;
   SetLength(Children._, NChildren);
   for i := 0 to NChildren - 1 do
      begin
      if UseCrossover then
         begin
         repeat
            IdParentA := RandAntibodyIndex(Parents, RankProp);
            IdParentB := RandAntibodyIndex(Parents, RankProp);
         until IdParentA <> IdParentB;
         Crossover(Children._[i], Parents._[IdParentA], Parents._[IdParentB]);
         end
      else
         begin
         IdParentA := RandAntibodyIndex(Parents, RankProp);
         AssignAntibody(Children._[i], Parents._[IdParentA]);
         end;
      Mutate(Children._[i], MutationRate);
      end;
   end;


// Replace a fraction of Antibodies with Children. 
// Which antibodies are replaced depends on the Replacement setting.
procedure Replace(
   var   Antibodies  :  TAntibodies;
   const Children    :  TAntibodies);
   var
         i, j, k,
         Which       :  Integer;
         Replaced    :  TBoolArray;
         Orders      :  array [0 .. 1] of TIntArray;
         Indices     :  TIntArray;
   begin
   k := Random(2);
   InitArray(Replaced, Antibodies.N, False);
   InitArray(Indices, {Len:} 2, {Value:} 0);
   Orders[0] := Antibodies.OrderFitness;
   Orders[1] := Antibodies.OrderSumW;
   for i := 0 to Children.N - 1 do
      begin
      case Replacement of
         repFitness, repSumWRankMul:
            j := Antibodies.OrderFitness[i];
         repDual:
            begin
            Which := Ord(((i + k) mod SumWFrac) = 1);
            repeat
               j := Orders[Which][Indices[Which]];
               Inc(Indices[Which]);
            until not Replaced[j];
            end;
         else
            Assert(False);
         end;
      AssignAntibody(Antibodies._[j], Children._[i]);
      Replaced[j] := True;
      //UpdateStats(Antibodies._[j], Antigens);
      end;
   end;


// Train Antigens on Antibodies for one generation using LearnRate
procedure OneGeneration(
   var   Antibodies  :  TAntibodies;
   const Antigens    :  TAntigens;
         LearnRate   :  Real);
   var
         Children    :  TAntibodies;
   begin
   Replicate(Children, Antibodies,
      {NChildren:} Round(LearnRate * Antibodies.N),
      {MutationRate:} 1 / (1 + Antigens.NVars), UseCrossover);
   Replace(Antibodies, Children);
   CalcStats(Antibodies, Antigens);
   end;


// Train Antigens on Antibodies for MaxGenerations
procedure Train(
   var   Antibodies     :  TAntibodies;
   const Antigens       :  TAntigens);
// const TrainParams    :  TTrainParams);
   var
         LearnRate      :  Real;
         i              :  Integer;
         BestAntibodies :  TAntibodies;
         FileStatus     :  Text;
         PathOut        :  AnsiString;
   const
         PathStatus     =  'Status.txt';
   begin
   // Initialize the antibodies
   BestAntibodies.Accuracy := 0;
   InitAntibodies(Antibodies, {N:} Round(Antigens.N * Abundance), Antigens);
   CalcStats(Antibodies, Antigens);
   
   OpenWrite(FileStatus, PathStatus);
   for i := 1 to MaxGens do
      begin
      // Train, update the best population
      LearnRate := LogBlend(1 / 2, 1 / Antibodies.N, (i - 1) / (MaxGens - 1));
      OneGeneration(Antibodies, Antigens, LearnRate);
      if Antibodies.Accuracy > BestAntibodies.Accuracy then
         CopyAntibodies(BestAntibodies, Antibodies);
         
      // Write the status, save antibodies for visualization
      WriteLn(FileStatus,
         i                    , Tab,
         LearnRate            , Tab,
         Antibodies.Coverage  , Tab,
         Antibodies.Accuracy  );
      if SaveVis then
         begin
         PathOut := PathVis + PadLeft(IntToStr(i), '0', 6) + '.txt';
         SaveAntibodiesVis(PathOut, Antibodies, Antigens);
         end;
      end;
   Close(FileStatus);
      
   CopyAntibodies(Antibodies, BestAntibodies);
   end;

{-----------------------<< Crossvalidation >>---------------------------------------}

// Return the Testset accuracy of a k-nearest classifier trained on Trainset
function kNNAccuracy(
   const Testset,
         Trainset    :  TAntigens;
         kNearest    :  Integer
         )           :  Real;
   var
         i, k, C,
         NCorrect    :  Integer;
         D           :  TRealArray;
         Votes       :  TRealArray;
         Order       :  TIntArray;
   begin
   SetLength(D, TrainSet.N);
   NCorrect := 0;
   for i := 0 to TestSet.N - 1 do
      begin
      // Pick the majority class of k nearest points
      SortByDistance(Order, Trainset, TestSet._[i].x);
      OrderRealArray(Order, D, soAScending);
      InitArray(Votes, TestSet.NClasses, 0);
      for k := 0 to kNearest - 1 do
         begin
         C := TrainSet._[Order[k]].Class_;
         Votes[C] := Votes[C] + 1;
         end;
      NCorrect := NCorrect + Ord(TestSet._[i].Class_ = RandMaxIndex(Votes));
      end;
   Result := NCorrect / TestSet.N;
   end;


// Run multiple rounds of k-fold crossvalidation, save results to FileResults.
// If kNN > 0, use kNN classification
procedure CrossValidate(
   const Antigens    :  TAntigens;
         Folds,
         Rounds      :  Integer;
   var   FileResults :  Text;
         kNN         :  Integer     =  0);
   var
         i, j, k,
         L, R        :  Integer;
         Order       :  TIntArray;
         TrainSet,
         TestSet     :  TAntigens;
         Antibodies  :  TAntibodies;
         Accuracies  :  TRealArrayN;
         FoldAccuracy,
         RoundCorrect,
         Mean, SD    :  Real;
   begin
   // Run multiple crossvalidation rounds with different data partitioning
   InitArrayN(Accuracies, 0);
   for j := 1 to Rounds do
      begin
      RandPerm(Order, Antigens.N, {Base:} 0);
      R := -1;
      RoundCorrect := 0;
      
      // Run k-fold crossvalidation
      for k := 1 to Folds do
         begin
         // Partition the antigens into the testing and the training set
         Write(k, ' ');
         L := R + 1;
         R := ((Antigens.N - 1) * k) div Folds;
         InitAntigens(TrainSet, Antigens.NVars, Antigens.NClasses);
         InitAntigens( TestSet, Antigens.NVars, Antigens.NClasses);
         for i := 0 to Antigens.N - 1 do
            begin
            if (i >= L) and (i <= R) then
               AppendAntigen(Testset,  Antigens._[Order[i]]) else
               AppendAntigen(Trainset, Antigens._[Order[i]]);
            end;
            
         // Get the classification accuracy
         if kNN = 0 then
            begin
            Train(Antibodies, TrainSet);
            FoldAccuracy := TestAccuracy(Antibodies, TestSet);
            end
         else
            FoldAccuracy := kNNAccuracy(Testset, Trainset, kNN);
         RoundCorrect := RoundCorrect + TestSet.N * FoldAccuracy;
         end;
         
      Append(Accuracies, RoundCorrect / Antigens.N);
      WriteLn;
      end;
      
   // Save the statistics
   SetTrueLength(Accuracies);
   GetMeanStandDev(Mean, SD, Accuracies._);
   WriteLn(FileResults, Mean, ' ', SD);
   end;


// Perform cross-validation on multiple datasets
procedure BatchCrossValidate;
   const
         NTests      =  6;
         Paths       :  array [1 .. NTests] of AnsiString
                     =  (  'Data_Iris.txt',
                           'Data_Wine.txt',
                           'Data_Glass.txt',
                           'Data_Crabs.txt',
                           'Data_Ionosphere.txt',
                           'Data_Sonar.txt'
                        );
         PathResults =  'Results.txt';
   var
         i           :  Integer;
         Antigens    :  TAntigens;
         FileResults :  Text;
   begin
   OpenWrite(FileResults, PathResults);
   for i := 1 to NTests do
      begin
      WriteLn(Paths[i]);
      LoadData(Antigens, DirData + PathDelim + Paths[i]);
      CrossValidate(Antigens, {Folds:} 5, {Rounds:} 6, FileResults, {kNN:} 0);
      end;
   Close(FileResults);
   end;

/////////////////////////////////////////////////////////////////////////////////////
var
      Antigens    :  TAntigens;
      Antibodies  :  TAntibodies;
begin
//Randomize;
LoadData(Antigens, DirData + PathDelim + 'Data_Iris.txt');
Train(Antibodies, Antigens);
SaveAntibodiesStats('Antibodies.txt', Antibodies);
//CrossValidate(Antigens, {Folds:} 5, {Rounds:} 6, {kNN:} 0);
//ReadLn;
//BatchCrossValidate;
end.

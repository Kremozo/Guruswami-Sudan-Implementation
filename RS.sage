class GSDecoder:
    def __init__(self, field, n, k, m=None):
        self.F = field
        self.n = n
        self.k = k
        self.R_xy = PolynomialRing(self.F, names=['x', 'y'])
        self.x, self.y = self.R_xy.gens()
        
        if m is None:
            self.m = 1
        else:
            self.m=m

        num_constraints = self.n * self.m * (self.m + 1) // 2
        D_approx = math.sqrt(2 * (self.k - 1) * num_constraints)
        self.l = math.floor(D_approx / (self.k - 1))
        
    def _get_degrees(self, n_constraints):
        monomials = []
        
        max_y_deg = max(self.l, self.m)
        max_x_deg = n_constraints 
        
        for j in range(max_y_deg+1):
            for i in range(max_x_deg+1):
                weight = i + (self.k - 1) * j
                monomials.append((weight, i, j))
        
        monomials.sort()
        
        return [ (m[1], m[2]) for m in monomials[:n_constraints + 1] ]

    def interpolate(self, points):
        """
        Find Q(x, y) passing through points with multiplicity m.
        """
        num_constraints = self.n * (self.m * (self.m + 1)) // 2
        
        powers = self._get_degrees(num_constraints)
        
        matrix_rows = []
        
        for (x_i, y_i) in points:
            # For each point, enforce that all derivatives of order < m vanish
            for r in range(self.m):
                for s in range(self.m - r):
                    row = []
                    for (i, j) in powers:
                        # If the monomial degree (i, j) is lower than the derivative order (r, s), the term vanishes (binomial is 0).
                        if i < r or j < s:
                            val = self.F(0)
                        else:
                            val = binomial(i, r) * binomial(j, s) * \
                                  (x_i**(i-r)) * (y_i**(j-s))
                        
                        row.append(val)
                    matrix_rows.append(row)
                    
        M = Matrix(self.F, matrix_rows)
        
        kernel = M.right_kernel()
        
        if kernel.dimension() == 0:
            return None
        
        coeffs = kernel.basis()[0]
        
        # Construct Q(x, y)
        Q = 0
        for idx, (i, j) in enumerate(powers):
            Q += coeffs[idx] * (self.x**i) * (self.y**j)
            
        return Q

    def factor_and_filter(self, Q):
        """
        Find roots y = f(x) of Q(x, y) such that deg(f) < k.
        """
        # View Q as a polynomial in y over F[x]
        R_x = PolynomialRing(self.F, 'x')
        R_y_over_x = PolynomialRing(R_x, 'y')
        
        try:
            # Convert multivariate Q to univariate y over x
            Q_uni = R_y_over_x(Q)
        except TypeError:
            print("Conversion to univariate failed.")
            return []

        # Find roots in F[x]
        roots = Q_uni.roots()
        
        valid_polynomials = []
        for (func, multiplicity) in roots:
            if func.degree() < self.k:
                valid_polynomials.append(func)
                
        return valid_polynomials

    def decode(self, received_points):
        Q = self.interpolate(received_points)
        if Q is None:
            print("Interpolation failed to find a polynomial.")
            return []
        
        #print(f"Interpolated Polynomial Q(x,y) found.")
        candidates = self.factor_and_filter(Q)
        return candidates


import time
import random
import math
import matplotlib.pyplot as plt

# --- PRELIMINARIES: ENCODING AND NOISE ---

def generate_message(k, F):
    """Generates a random message polynomial of degree < k."""
    R = PolynomialRing(F, 'x')
    return R.random_element(degree=k-1)

def encode_rs(message_poly, points, n):
    """Encodes the message by evaluating it at n distinct points."""
    codeword = []
    for i in range(n):
        val = message_poly(points[i])
        codeword.append((points[i], val))
    return codeword

def add_errors(codeword, num_errors, F):
    """Corrupts 'num_errors' y-values in the codeword."""
    n = len(codeword)
    corrupted = list(codeword) # Copy
    
    error_indices = random.sample(range(n), num_errors)
    
    for idx in error_indices:
        x_val, original_y = corrupted[idx]
        
        while True:
            noise = F.random_element()
            if noise != original_y:
                corrupted[idx] = (x_val, noise)
                break
                
    return corrupted

# --- EXPERIMENT SUITE CLASS ---

class GSExperimentSuite:
    def __init__(self, n, k, field_size=256):
        # 1. Setup Field and Code Parameters
        self.F = GF(field_size, 'a')
        self.points = self.F.list()[:n] # Evaluation points
        self.n = n
        self.k = k
        self.d_min = n - k + 1
        
        # Calculate Bounds
        self.unique_bound = (n - k) // 2
        self.johnson_bound = n - int(math.sqrt(n * (k - 1)))
        
        print(f"--- Code Parameters RS({n}, {k}) over GF({field_size}) ---")
        print(f"Unique Decoding Radius: {self.unique_bound}")
        print(f"Johnson Radius: {self.johnson_bound}")
        print("-------------------------------------------------------")

    def run_radius_experiment(self, ms=[1, 2, 3], trials=20):
        """
        Experiment 1: Decoding Radius Analysis
        - varies multiplicity (m)
        - varies error count (e)
        - Logs detailed results to 'experiment_logs.txt' and console.
        - Generates 'exp1_radius_comparison.png'
        """
        filename = "experiment_logs.txt"
        
        with open(filename, "w") as f:
            
            def log(text):
                print(text)
                f.write(text + "\n")

            log(f"{'='*60}")
            log(f" EXPERIMENT 1: DECODING RADIUS ANALYSIS")
            log(f" Parameters: N={self.n}, K={self.k}, Field={self.F}")
            log(f" Multiplicities tested: {ms}")
            log(f"{'='*60}\n")
            
            start_e = int(self.unique_bound - 2)
            end_e = int(self.johnson_bound + 3)
            error_counts = range(start_e, end_e)
            
            plt.figure(figsize=(10, 6))
            
            for m_val in ms:
                log(f"\n{'#'*60}")
                log(f" STARTING BATCH: Multiplicity m={m_val}")
                log(f"{'#'*60}")
                
                success_rates = []
                decoder = GSDecoder(self.F, self.n, self.k, m=m_val)
                for e in error_counts:
                    log(f"\n[TEST] Testing Error Count e={e}...")
                    successes = 0
                    
                    for t in range(1, trials + 1):
                        msg = generate_message(self.k, self.F)
                        codeword = encode_rs(msg, self.points, self.n)
                        received = add_errors(codeword, e, self.F)
                        
                        candidates = decoder.decode(received)
                        
                        found = msg in candidates
                        if found:
                            successes += 1
                            status = "FOUND "
                        else:
                            status = "FAILED"
                        
                        log(f"   Trial {t:02d}: {status} (Errors: {e})")
                        log(f"             Input Poly:   {msg}")
                        
                        if len(candidates) > 0:
                            cand_str = str(candidates).replace('\n', ' ') 
                            log(f"             Output List:  {cand_str}")
                        else:
                            log(f"             Output List:  [ EMPTY ]")
                        
                        log("-" * 40)

                    rate = float(successes / trials)
                    success_rates.append(rate)
                    log(f">> SUMMARY for e={e}: {int(rate*100)}% Success ({successes}/{trials})")

                plt.plot(list(error_counts), success_rates, marker='o', label=f'm={m_val}')

            plt.axvline(x=float(self.unique_bound), color='r', linestyle='--', label='Unique Bound')
            plt.axvline(x=float(self.johnson_bound), color='g', linestyle='--', label='Johnson Bound')
            
            plt.title(f"Decoding Radius vs Multiplicity (m)\nRS({self.n}, {self.k}) over GF(2^8)")
            plt.xlabel("Number of Errors")
            plt.ylabel("Success Probability")
            plt.legend()
            plt.grid(True, which='both', linestyle='--', alpha=0.7)
            
            plot_filename = "exp1_radius_comparison.png"
            plt.savefig(plot_filename)
            log(f"\n[INFO] Plot saved to {plot_filename}")
            log(f"[INFO] Full logs saved to {filename}")

    def run_runtime_experiment(self, max_m=5, trials=10):
        """
        Exp 2: Runtime vs Multiplicity (m)
        Demonstrates the cost of the 'heavy linalg'.
        """
        print(f"\n[Experiment 2] Running Runtime Analysis...")
    
        errors = self.johnson_bound 
        
        m_values = range(1, max_m + 1)
        avg_times = []
        
        for m in m_values:
            decoder = GSDecoder(self.F, self.n, self.k, m=m)
            total_time = 0
            
            for _ in range(trials):
                msg = generate_message(self.k, self.F)
                codeword = encode_rs(msg, self.points, self.n)
                received = add_errors(codeword, errors, self.F)
                
                start = time.time()
                decoder.decode(received)
                total_time += (time.time() - start)
            
            avg = total_time / trials
            avg_times.append(avg)
            print(f"  Multiplicity m={m} | Avg Time: {avg:.4f}s")

        # Plotting
        plt.figure(figsize=(8, 5))
        plt.plot(m_values, avg_times, marker='s', color='orange')
        plt.title("Decoding Runtime vs Multiplicity Parameter")
        plt.xlabel("Multiplicity (m)")
        plt.ylabel("Time (seconds)")
        plt.grid(True)
        plt.savefig("exp2_runtime.png")
        print("-> Saved plot to exp2_runtime.png")

    def run_list_size_experiment(self, m=3, trials=100):
        """
        Exp 3: List Size Distribution
        How many 'ghost' polynomials do we actually find?
        """
        print(f"\n[Experiment 3] Analyzing List Sizes...")
        errors = 29
        list_sizes = []
        
        decoder = GSDecoder(self.F, self.n, self.k, m=m)
        
        for _ in range(trials):
            msg = generate_message(self.k, self.F)
            codeword = encode_rs(msg, self.points, self.n)
            received = add_errors(codeword, errors, self.F)
            
            candidates = decoder.decode(received)
            list_sizes.append(len(candidates))

        plt.figure(figsize=(8, 5))
        plt.hist(list_sizes, bins=range(min(list_sizes), max(list_sizes) + 2), align='left', rwidth=0.8)
        plt.title(f"Histogram of List Sizes (Errors={errors})")
        plt.xlabel("List Size")
        plt.ylabel("Frequency")
        plt.xticks(range(min(list_sizes), max(list_sizes) + 1))
        plt.savefig("exp3_listsize.png")
        print("-> Saved plot to exp3_listsize.png")

    def run_rate_stress_test(self, trials=20):
        """
        Exp 4 (Enhanced): High Rate vs Low Rate Analysis.
        Tests if we can decode BEYOND the Unique bound across different K values.
        Run multiple trials to get statistical confidence and track list sizes.
        """
        print(f"\n[Experiment 4] Running Rate Stress Test (Pathological K)...")
        print(f"Testing at ~50% of the way to Johnson Bound with m=2 over {trials} trials.")
        
        k_values = [4, 16, 32, 60] 
        
        print(f"{'K':<5} | {'Unique':<7} | {'Johnson':<7} | {'Target':<7} | {'Success%':<9} | {'Avg |L|':<8}")
        print("-" * 65)
        
        for k_curr in k_values:
            unique = (self.n - k_curr) // 2
            johnson = self.n - int(math.sqrt(self.n * (k_curr - 1)))
            gain = johnson - unique
            
            if gain < 2:
                errors_to_test = unique
            else:
                errors_to_test = unique + (gain // 2)
            
            decoder = GSDecoder(self.F, self.n, k_curr, m=2)
            
            success_count = 0
            total_list_size = 0
            
            for _ in range(trials):
                msg = generate_message(k_curr, self.F)
                codeword = encode_rs(msg, self.points, self.n)
                received = add_errors(codeword, errors_to_test, self.F)
                
                candidates = decoder.decode(received)
                
                if msg in candidates:
                    success_count += 1
                
                total_list_size += len(candidates)
            
            success_rate = float((success_count / trials) * 100)
            avg_list_size = float(total_list_size / trials)
            
            print(f"{k_curr:<5} | {unique:<7} | {johnson:<7} | {errors_to_test:<7} | {success_rate:<8.1f}% | {avg_list_size:<8.2f}")

    def run_interleaved_collision(self, trials=20):
        """
        Exp 5 (Randomized): 'Signal Collision' / Random Interleaving.
        Simulates two codewords 'colliding' on the channel.
        For each position, we randomly pick the symbol from C1 or C2.
        """
        print(f"\n[Experiment 5] Running Random Interleaved Collision (Signal Jamming)...")
        
        k_col = 4
        m_col = 5 
        
        decoder = GSDecoder(self.F, self.n, k_col, m=m_col)
        
        success_both = 0
        avg_list_size = 0
        
        print(f"{'Trial':<6} | {'Msg1 Found':<10} | {'Msg2 Found':<10} | {'List Size'}")
        print("-" * 50)
        
        for t in range(trials):
            msg1 = generate_message(k_col, self.F)
            msg2 = generate_message(k_col, self.F)
            while msg1 == msg2:
                msg2 = generate_message(k_col, self.F)
            
            c1 = encode_rs(msg1, self.points, self.n)
            c2 = encode_rs(msg2, self.points, self.n)
            
            received = []
            for i in range(self.n):
                # Flip a coin: 50% chance to take from C1, 50% from C2
                if random.random() < 0.5:
                    received.append(c1[i])
                else:
                    received.append(c2[i])
            
            candidates = decoder.decode(received)
            
            found1 = msg1 in candidates
            found2 = msg2 in candidates
            l_size = len(candidates)
            
            avg_list_size += l_size
            if found1 and found2:
                success_both += 1
                
            if t < 5:
                print(f"{t+1:<6} | {str(found1):<10} | {str(found2):<10} | {l_size}")
        
        final_avg_list = float(avg_list_size / trials)
        double_capture_rate = float((success_both / trials) * 100)

        print("-" * 50)
        print(f"Total Trials: {trials}")
        print(f"Double Capture Rate: {double_capture_rate:.1f}%")
        print(f"Average List Size: {final_avg_list:.2f}")
# --- EXECUTION BLOCK ---

# Parameters: RS(64, 16) over GF(256)
# This gives a unique bound of (64-16)/2 = 24 errors
# Johnson bound is approx 64 - sqrt(64*15) ~ 64 - 30.9 = 33 errors
if __name__ == "__main__":
    suite = GSExperimentSuite(n=64, k=16)
    
    # Run the experiments
    suite.run_radius_experiment(ms=[1, 2, 3,4], trials=20)
    suite.run_runtime_experiment(max_m=5, trials=20)
    suite.run_list_size_experiment(m=3, trials=100)
    suite.run_rate_stress_test()
    suite.run_interleaved_collision(trials=20)

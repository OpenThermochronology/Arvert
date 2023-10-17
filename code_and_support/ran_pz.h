struct Ranq1 {
	Ullong v;
	Ranq1(Ullong j) : v(4101842887655102017LL) {
		v ^= j;
		v = int64();
	}
	inline Ullong int64() {
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v * 2685821657736338717LL;
	}
	inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};